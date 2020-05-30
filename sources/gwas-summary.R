#!/usr/bin/Rscript

# INFO   : Create different summarized reports (tables, plots, Venn diagrams) for multiGWAS tool analysis
# AUTHOR : Luis Garreta (lgarreta@agrosavia.co)
# DATA   : feb/2020
# LOG: 
#	r5.0: Heuristic to select best gene action model
#	r4.1: Fixed transparency and axis error (options (bitmapType="cairo"))
#	r4.0: Modified to work with markdown, but better to only report outputs (PNGs) to be included by markdown
#	r3.0: Manhattan and QQ plots. Formated to create a Markdown report (not yet)
#	r2.0: Improve selection of data from tables. Titles to graphics and files
#	r1.0: Message to log files
#	r0.9: Full working with funtions: create snpTables and ven diagrams using parameters
#	r0.8: Create venn diagrams, summary table of first Ns"
#     


suppressMessages (library (dplyr))
suppressMessages (library (qqman))
suppressMessages (library (VennDiagram))
suppressMessages (library (config))  # For read config file


#source ("gwas-heatmap.R")            # Module with functions to create summaries: tables and venn diagrams
#source ("gwas-preprocessing.R")      # Module with functions to convert between different genotype formats 

options (bitmapType="cairo", width=300)
#options(scipen=999)

#-------------------------------------------------------------
# Main function
# Input files are taken from input dir
# Outupt are written to output dir
#-------------------------------------------------------------
main <- function () {

	msgmsg ("Main...")

	inputDir     = "out/"
	genotypeFile  = "out/filtered-gwasp4-genotype.tbl"
	phenotypeFile = "out/filtered-gwasp4-phenotype.tbl"
	outputDir   = "report/"
	gwasModel    = "Full"
	nBest = 10

	createReports (inputDir, genotypeFile, phenotypeFile, 
				   gwasModel, outputDir, nBest)
}

#-------------------------------------------------------------
# Function to create different reports:
#	1- 1 table of best SNPs
#	2- 1 table of significative SNPs
#   3- 1 Venn diagram of best SNPs
#   4- 1 Venn diagram of significative SNPs
#	5- 1 multiplot of 4x4 manhattan and QQ plots
#-------------------------------------------------------------
# Only for Manfred's data
#-------------------------------------------------------------
reduceSNPNames <- function (resultFiles) {
	for (filename in resultFiles) {
		data = read.csv (file=filename, sep="\t")
		if (grepl ("TASSEL", filename) || grepl ("GWASpoly", filename))
			data$Marker = gsub ("solcap_snp_","", data$Marker)
		else
			data [,"SNP"] = gsub ("solcap_snp_","", data$SNP)
		#newFilename = paste0 (strsplit (filename,"[.]")[[1]][1], "-new.csv")
		write.table (file=filename, data, sep="\t", quote=F, row.names=F)
	}
}
#-------------------------------------------------------------
#-------------------------------------------------------------
createReports <- function (inputDir, genotypeFile, phenotypeFile, gwasModel, outputDir, nBest) 
{
	message ("Starting summary...")
	# Define filenams for outputs
	fileBestScores               = paste0(outputDir,  "/out-multiGWAS-scoresTable-best.scores")
	fileSignificativeScores      = paste0(outputDir,  "/out-multiGWAS-scoresTable-significatives.scores")
	fileBestVennDiagram          = paste0 (outputDir, "/out-multiGWAS-vennDiagram-best")
	fileSignificativeVennDiagram = paste0 (outputDir, "/out-multiGWAS-vennDiagram-significatives")
	fileManhattanPlotPNG         = paste0 (outputDir, "/out-multiGWAS-manhattanQQ-plots.png")
	fileManhattanPlotPDF         = paste0 (outputDir, "/out-multiGWAS-manhattanQQ-plots.pdf")


	msgmsg ("Creating reports for ", gwasModel, " model...")
	createDir (outputDir)

	msgmsg ("Writing table input config parameters...")
	writeConfigurationParameters (inputDir, outputDir)

	# Get filenames of results for each of the four GWAS tools
	resultFiles =  list.files(inputDir, pattern=sprintf ("(^(tool).*(%s).*[.](csv))", gwasModel), full.names=T)
	reduceSNPNames (resultFiles)

	msgmsg ("Creating table with summary results...")
	snpTables = markersSummaryTable (resultFiles, gwasModel, nBest)

	msgmsg ("Writing table with ", nBest, " best ranked SNPs Table...")
	write.table (file=fileBestScores, snpTables$best, row.names=F,quote=F, sep="\t")

	msgmsg ("Writing table with significative SNPs...")
	write.table (file=fileSignificativeScores, snpTables$significatives, row.names=F,quote=F, sep="\t")

	msgmsg ("Writing Venn diagram with best SNPs...")
	commonBest = markersVennDiagrams (snpTables$best, gwasModel, "Best", fileBestVennDiagram)

	msgmsg ("Writing Venn diagram with significative SNPs...")
	commonSign = markersVennDiagrams (snpTables$significatives, gwasModel, "Significatives", fileSignificativeVennDiagram)

	msgmsg ("Writing Manhattan and QQ plots...")
	png (fileManhattanPlotPNG, width=11, height=15, units="in", res=120)
	op=markersManhattanPlots (resultFiles, gwasModel, commonBest, commonSign, snpTables, outputDir, nBest)
	dev.off()

	pdf (fileManhattanPlotPDF, width=11, height=15)
	op=markersManhattanPlots (resultFiles, gwasModel, commonBest, commonSign, snpTables, outputDir, nBest)
	par (op)
	dev.off()

	
	# Create heat maps
	msgmsg ("Creating SNP heatmaps for the best ranked SNPs...")
	genoNumericFilename = ACGTToNumericGenotypeFormat (genotypeFile)
	createHeatmapForSNPList (outputDir, genotypeFile, genoNumericFilename, phenotypeFile, commonBest)
}
#-------------------------------------------------------------
# Calculate the inflation factor from -log10 values
# It can fire warning, here they are hidign
#-------------------------------------------------------------
calculateInflationFactor <- function (scores)
{
	oldw <- getOption("warn")
	options(warn = -1)

	remove <- which(is.na(scores))
	if (length(remove)>0) 
		x <- sort(scores[-remove],decreasing=TRUE)
	else 
		x <- sort(scores,decreasing=TRUE)

	pvalues = 10^-x
	chisq <- na.omit (qchisq(1-pvalues,1))
	delta  = round (median(chisq)/qchisq(0.5,1), 3)

	options (warn = oldw)

	return (list(delta=delta, scores=x))
}

#------------------------------------------------------------------------
#------------------------------------------------------------------------
markersManhattanPlots <- function (resultFiles, gwasModel, commonBest, commonSign, snpTables, outputDir, nBest) {
	#resultFiles =  list.files(inputDir, pattern=sprintf ("(^(tool).*(%s).*[.](csv))", gwasModel), full.names=T)
	op <- par(mfrow = c(4,2), mar=c(3.5,3.5,3,1), oma=c(0,0,0,0), mgp = c(2.2,1,0))
	for (filename in resultFiles) {
		data           = read.table (file=filename, header=T)
		data           = data [!is.na (data$P),]

		names          = unlist (strsplit (basename (filename), "[-|.]"))
		mainTitle      = paste0 (names[2],"-", names [3])

		if (grepl ("GWASpoly", filename)) {
			tool = "GWASpoly"
			data = selectBestModel (data, nBest, tool)
			gwasResults = data.frame (SNP=data$Marker, CHR=data$Chrom, BP=data$Position, P=10^-data$SCORE)
		}
		else if (grepl ("SHEsis", filename)) {
			tool = "SHEsis"
			gwasResults = data.frame (SNP=data$SNP, CHR=data$CHR, BP=data$POS, P=10^-data$SCORE)
		}
		else if (grepl ("PLINK", filename)) {
			tool = "PLINK"
			gwasResults = data.frame (SNP=data$SNP, CHR=data$CHR, BP=data$POS, P=10^-data$SCORE)
		}
		else if (grepl ("TASSEL", filename)) {
			tool = "TASSEL"
			data = selectBestModel (data, nBest, tool)
			gwasResults = data.frame (SNP=data$Marker, CHR=data$Chr, BP=data$Pos, P=10^-data$SCORE)
		}

		# Select limiting lines for Manhattan (best and significant)
		bestThresholdScore = data [nBest,c("SCORE")]
		bestThreshold      = 10^-bestThresholdScore
		signThresholdScore = data [1, "THRESHOLD"]
		signThreshold      = 10^-signThresholdScore

		# 
		ss = snpTables$significatives
		if (tool %in% ss$TOOL)
			signThresholdScore = min (ss [ss$TOOL==tool,"SCORE"])
		else
			signThresholdScore = ceiling (data[1, "SCORE"]) 

		bestSNPsTool     = unlist (select (filter (snpTables$best, TOOL==tool), "SNP"))
		sharedSNPs       = intersect (commonBest, bestSNPsTool)
		colorsBlueOrange = c("blue4", "orange3")
		ylims   = c (0, ceiling (signThresholdScore))
		manhattan(gwasResults,col = c("orange", "midnightblue"), highlight=sharedSNPs, annotatePval=bestThreshold, annotateTop=F,
				  suggestiveline=bestThresholdScore, genomewideline=signThresholdScore, main=mainTitle, logp=T, cex=2)

		text (x=0, y=signThresholdScore*0.92, "               Significants",, col="red", pos=4)
		text (x=0, y=bestThresholdScore*0.92, "Best",, col="blue", pos=4)

		datax = calculateInflationFactor (-log10 (gwasResults$P))
		qq (gwasResults$P)
		mtext (bquote(lambda[GC] == .(datax$delta)), side=3, line=-2, cex=0.7)
		#title (bquote(lambda[GC] == .(datax$delta)))
	}
	#par (op)
	#dev.off()
	return (op)
}

#-----------------------------------------------------------
# Select best N SNPs from multiple action models (for GWASpoly and TASSEL)
# Main criteria is GC
# PLINK also can produce info of more action models using options
#-----------------------------------------------------------
selectBestModel <- function (data, nBest, tool) {
	outFilename = paste0 ("report/model-scores-", tool, ".csv")
	nBest = 50
	# Select main columns
	dr = data [,c("Marker","GC","MODEL","SCORE", "THRESHOLD", "DIFF")]; 

	# Order by nBest, DIFF, GC
	do = dr [order (dr$MODEL,-dr$DIFF),]; 

	# Reduce to groups of nBest
	dm = Reduce (rbind, by(do, do["MODEL"], head, n=nBest)); 

	# Add Count of SNPs
	summ   = data.frame (add_count (dm, Marker, sort=T)); 
	summMdl = aggregate (x=summ$n, by=list(MODEL=summ$MODEL, GC=summ$GC), 	FUN=sum)

	# Calculate best model score
	totalNs    = length (summMdl$MODEL) * nBest
	modelScore = summMdl$x/totalNs  + 1 - abs (1-summMdl$GC)
	summScores = cbind (summMdl, score=modelScore)
	summScores = summScores [order (summScores$score, summScores$MODEL, decreasing=T),]
	write.table (file=outFilename, summScores, quote=F, row.names=F, sep="\t")
	bestModel = summScores [1, "MODEL"]

	#message ("Best model is: ", bestModel)

	# Select SNPs for model and sort by DIFF
	dataModel = data [data[,"MODEL"] %in% bestModel,]
	dataModel = dataModel [order (-dataModel$DIFF),]

	return (dataModel)
}

#------------------------------------------------------------------------
# Create Venn diagram of common markers using info from summary table
#------------------------------------------------------------------------
markersVennDiagrams <- function (summaryTable, gwasModel, scoresType, outFile){
	flog.threshold(ERROR)
	# Set lists for Venn diagrams
	x <- list()
	x$GWASpoly = summaryTable %>% filter (TOOL %in% "GWASpoly") %>% select (SNP) %>% .$SNP
	x$SHEsis   = summaryTable %>% filter (TOOL %in% "SHEsis") %>% select (SNP) %>% .$SNP
	x$PLINK    = summaryTable %>% filter (TOOL %in% "PLINK")  %>% select (SNP) %>% .$SNP
	x$TASSEL   = summaryTable %>% filter (TOOL %in% "TASSEL") %>% select (SNP) %>% .$SNP

	# Create Venn diagram
	mainTitle = paste0(gwasModel, "-", scoresType)
	COLORS= c("red", "blue", "yellow", "green")
	v0 <- venn.diagram(x, height=3000, width=3000, alpha = 0.5, filename = NULL, # main=mainTitle,
						col = COLORS, cex=0.9, margin=0.0,
						fill = COLORS)

	overlaps <- calculate.overlap(x)
	overlaps <- rev(overlaps)

	posOverlap = as.numeric (gsub ("a","", (names (overlaps))))
	for (i in 1:length(overlaps)){
		pos = posOverlap [i]
		v0[[pos+8]]$label <- paste(overlaps[[i]], collapse = "\n")
 	}

	WIDTH  = 9
	HEIGHT = 12

 	png (paste0 (outFile,".png"), width=WIDTH, height=HEIGHT, units="in", res=120)
	grid.draw(v0)
	dev.off()
	pdf (paste0 (outFile,".pdf"), width=WIDTH,height=HEIGHT)
	grid.draw(v0)
	dev.off()
	
	# Get shared SNPs
	dataSNPsNs     = data.frame (add_count (summaryTable, SNP, sort=T)); 
	dataSNPsShared = dataSNPsNs[dataSNPsNs$n > 1,]
	dataSNPsNoDups = dataSNPsShared [!duplicated (dataSNPsShared$SNP),]
	sharedSNPs     = dataSNPsNoDups$SNP


	return (sharedSNPs)
}

#------------------------------------------------------------------------
# Create a summary table of best and significative markers
#------------------------------------------------------------------------
markersSummaryTable <- function (resultFiles, gwasModel, nBest) {
	#files =  list.files(inputDir, pattern=sprintf ("(^(tool).*(%s).*[.](csv))", gwasModel), full.names=T)
	msgmsg ("Creating summary table...")
	summaryTable = data.frame ()

	tool=""
	for (f in resultFiles) {
		data <- read.table (file=f, header=T)
		#if (nrow(data)>nBest) data=data [1:nBest,] 

		flagNewData = F
		if (grepl("GWASpoly", f)) {
			tool    = "GWASpoly"
			data    = selectBestModel (data, nBest, tool)
			chrom   = data$Chrom
			pos	    = data$Position
			snps    = data$Marker
			flagNewData = T
		}else if (grepl ("PLINK", f)) {
			tool    = "PLINK"
			chrom   = data$CHR
			pos	    = data$POS
			snps    = data$SNP
			flagNewData = T
		}else if (grepl ("TASSEL", f)) {
			tool    = "TASSEL"
			data    = selectBestModel (data, nBest, tool)
			chrom   = data$Chr
			pos		= data$Pos
			snps    = data$Marker
			flagNewData = T
		}else if (grepl ("SHEsis", f)) {
			tool    = "SHEsis"
			chrom   = data$CHR
			pos     = data$POS
			snps    = data$SNP
			flagNewData = T
		}

		# Set values with general column names
		model   = data$MODEL
		gcs     = data$GC
		pVal	= data$P
		pscores = data$SCORE
		tscores = data$THRESHOLD
		signf   = pscores >= tscores

		if (flagNewData==T) {
			dfm = data.frame (TOOL=tool, MODEL=model, GC=gcs, SNP=snps, CHROM=chrom, POSITION=pos, 
							  PVALUE = round (pVal,6), SCORE=round (pscores, 4), THRESHOLD=round (tscores,4), SIGNIFICANCE=signf )
			#dfm = dfm %>% distinct (SNP, .keep_all=T)
			dfm = dfm [!duplicated (dfm$SNP),]
			if (nrow(dfm)>nBest) dfm=dfm [1:nBest,] 
			summaryTable <- rbind (summaryTable, dfm)
			flagNewData = F
		}
	}

	summaryTable = summaryTable [which(!is.na(summaryTable$SIGNIFICANCE)),]
	summarySignificatives = summaryTable %>% filter (SIGNIFICANCE%in%T) 

	return (list (best=summaryTable, significatives=summarySignificatives))
}

#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
hd <- function (data, m=10,n=10, tool="") {
	filename = paste0 ("x",tool,"-", deparse (substitute (data)),".csv")
	msgmsg (deparse (substitute (data)),": ", dim (data))
	if (is.null (dim (data)))
		print (data [1:10])
	else if (ncol (data) < 10) 
		print (data[1:m,])
	else if (nrow (data) < 10)
		print (data[,1:n])
	else 
		print (data [1:m, 1:n])

	write.table (file=filename, data, quote=F, sep="\t", row.names=F)
}
#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msgmsg <- function (...) 
{
		messages = unlist (list (...))
		cat (">>>>", messages, "\n")
}

#-------------------------------------------------------------
# Add label to filename
#-------------------------------------------------------------
addLabel <- function (filename, label, newExt=NULL)  {
	nameext = strsplit (filename, split="[.]")
	name    = nameext [[1]][1] 
	if (is.null (newExt))
		ext     = nameext [[1]][2] 
	else
		ext     = newExt
	newName = paste0 (nameext [[1]][1], "-", label, ".", ext )
	return (newName)
}

#----------------------------------------------------------
# Create dir, if it exists the it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		name  = basename (newDir)
		path  = dirname  (newDir)
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", path, name)
			if (dir.exists (oldDir) == T) {
				checkOldDir (oldDir)
			}

			file.rename (newDir, oldDir)
		}
	}

	checkOldDir (newDir)
	system (sprintf ("mkdir %s", newDir))
}

#-------------------------------------------------------------
# Get params from config file and define models according to ploidy
#-------------------------------------------------------------
## @knitr writeConfigurationParameters
writeConfigurationParameters <- function (inputDir, outputDir) 
{
	configFile = paste0(inputDir, list.files (inputDir, pattern="config")[1])

	params = config::get (file=configFile) 

	paramsDF = data.frame (PARAMETER=character(), VALUE=character ())
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Genotype filename", VALUE=toString (params$genotypeFile)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Phenotype filename", VALUE=toString (params$phenotypeFile)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Significance level (Genome-wide significance level)", VALUE=toString (params$significanceLevel)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Correction method (Bonferroni or FDR)", VALUE=toString (params$correctionMethod)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="GWAS model (Full or Naive)", VALUE=toString (params$gwasModel)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="nBest (Number of best-ranked SNPs to be reported)", VALUE=toString (params$nBest)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="Filtering (TRUE or FALSE)", VALUE=toString (params$filtering)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="MIND Filter (Individual with missing genotype)", VALUE=toString (params$MIND)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="GENO Filter (SNPs with missing genotype)", VALUE=toString (params$GENO)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="MAF Filter (Minor allele frequency)", VALUE=toString (params$MAF)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="HWE Filter (Hardy-Weinberg test)", VALUE=toString (params$HWE)))
	paramsDF = rbind  (paramsDF, data.frame (PARAMETER="GWAS Tools", VALUE=toString (params$tools)))

	outName = paste0(outputDir, "/out-multiGWAS-inputParameters.tbl")
	write.table (file=outName, paramsDF, quote=F, sep="\t", row.names=F)
	return (paramsDF)
}
#-------------------------------------------------------------
# Get alternate allele from a row of alleles
#-------------------------------------------------------------
#-------------------------------------------------------------
# Call to main function (first lines)
#-------------------------------------------------------------
#main ()

