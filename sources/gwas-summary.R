#!/usr/bin/Rscript

# INFO   : Create different summarized reports (tables, plots, Venn diagrams) for multiGWAS tool analysis
# AUTHOR : Luis Garreta (lgarreta@agrosavia.co)
# DATA   : feb/2020
# LOG: 
#	r5.1: Fixed chord diagram when no shared SNPs
#	r5.0: Heuristic to select best gene action model
#	r4.1: Fixed transparency and axis error (options (bitmapType="cairo"))

HOME = Sys.getenv ("MULTIGWAS_HOME")
.libPaths (paste0(HOME, "/opt/Rlibs"))

suppressMessages (library (dplyr))
suppressMessages (library (qqman))
suppressMessages (library (VennDiagram))
suppressMessages (library (config))          # For read config file
suppressMessages (library ("RColorBrewer"))  # For chord diagrams
suppressMessages (library(circlize))         # For chord diagrams


options (bitmapType="cairo", width=300)
#options(scipen=999)

#-------------------------------------------------------------
# Main function
# Input files are taken from input dir
# Outupt are written to output dir
#-------------------------------------------------------------
main <- function () {
	options (warn=1)
	#source ("lglib01.R")
	source ("gwas-heatmap.R")            # Module with functions to create summaries: tables and venn diagrams
	source ("gwas-preprocessing.R")      # Module with functions to convert between different genotype formats 

	msgmsg ("Main...")

	inputDir     = "out/"
	genotypeFile  = "out/filtered-gwasp4-genotype.tbl"
	phenotypeFile = "out/filtered-gwasp4-phenotype.tbl"
	outputDir   = "report/"
	gwasModel    = "Naive"
	nBest  = 5
	ploidy = 4

	createReports (inputDir, genotypeFile, phenotypeFile, 
				   ploidy, gwasModel, outputDir, nBest)
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
reduceSNPNames <- function (listOfResultsFile) {
	for (filename in listOfResultsFile) {
		data = read.csv (file=filename, sep="\t", check.names=F)
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
createReports <- function (inputDir, genotypeFile, phenotypeFile, ploidy, gwasModel, outputDir, nBest) 
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
	config = writeConfigurationParameters (inputDir, outputDir)

	# Get filenames of results for each of the four GWAS tools
	listOfResultsFile =  list.files(inputDir, pattern=sprintf ("(^(tool).*(%s).*[.](csv))", gwasModel), full.names=T)
	reduceSNPNames (listOfResultsFile)

	msgmsg ("Creating table with summary results...")
	snpTables = markersSummaryTable (listOfResultsFile, gwasModel, nBest)

	msgmsg ("Writing table with ", nBest, " best ranked SNPs Table...")
	write.table (file=fileBestScores, snpTables$best, row.names=F,quote=F, sep="\t")

	msgmsg ("Writing table with significative SNPs...")
	write.table (file=fileSignificativeScores, snpTables$significatives, row.names=F,quote=F, sep="\t")

	msgmsg ("Writing Venn diagram with best SNPs...")
	commonBest = markersVennDiagrams (snpTables$best, gwasModel, "Best", fileBestVennDiagram)

	msgmsg ("Writing Venn diagram with significative SNPs...")
	commonSign = markersVennDiagrams (snpTables$significatives, gwasModel, "Significatives", fileSignificativeVennDiagram)

	msgmsg ("Writing Manhattan and QQ plots...")
	png (fileManhattanPlotPNG, width=11, height=15, units="in", res=90)
	op=markersManhattanPlots (listOfResultsFile, gwasModel, commonBest, commonSign, snpTables, outputDir, nBest)
	dev.off()

	pdf (fileManhattanPlotPDF, width=11, height=15)
	op=markersManhattanPlots (listOfResultsFile, gwasModel, commonBest, commonSign, snpTables, outputDir, nBest)
	par (op)
	dev.off()


	# Create heat maps
	msgmsg ("Creating SNP heatmaps for the best ranked SNPs...")
	genoNumericFilename = ACGTToNumericGenotypeFormat (genotypeFile, ploidy)

	createHeatmapForSNPList (outputDir, genotypeFile, genoNumericFilename, phenotypeFile, commonBest)

	# Create chord diagrams
	msgmsg ("Creating chord diagrams for chromosome vs SNPs...")
	createChordDiagramSharedSNPs (fileBestScores)

	# Call to rmarkdown report
	createMarkdownReport (config)

}

#-------------------------------------------------------------
# Create Rmkardown report
#-------------------------------------------------------------
createMarkdownReport  <- function (config) {
	msg ("Creating html rmarkdown report...")
	outputFile = paste0 (config$workingDir, "/multiGWAS-report.html") 
	title      = paste0 ("MultiGWAS report for ", config$gwasModel, " GWAS model")


	# Create html with embbeded images (using javascrips) for Java WebView
	rmarkdown::render (paste0(HOME,"/sources/gwas-markdown.Rmd"), output_file=outputFile, output_format="html_document", 
					   params=list (workingDir=config$workingDir, reportTitle=title, nBest=config$nBest), quiet=T)

	# Create html with external images (using subfolders) for Java JEditorPane
	#rmarkdown::render (paste0(HOME,"/sources/gwas-markdown.Rmd"), output_file=outputFile, 
	#				   output_format="html_document", output_options=list(self_contained=F),
	#				   params=list (workingDir=config$workingDir, reportTitle=title, nBest=config$nBest), quiet=T)
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
markersManhattanPlots <- function (listOfResultsFile, gwasModel, commonBest, commonSign, snpTables, outputDir, nBest) {
	#listOfResultsFile =  list.files(inputDir, pattern=sprintf ("(^(tool).*(%s).*[.](csv))", gwasModel), full.names=T)
	op <- par(mfrow = c(4,2), mar=c(3.5,3.5,3,1), oma=c(0,0,0,0), mgp = c(2.2,1,0))
	for (filename in listOfResultsFile) {
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
		#title (bquote(lambda[GC] == .(data$delta)))
	}
	#par (op)
	#dev.off()
	return (op)
}

#-----------------------------------------------------------
# Select best N SNPs from multiple action models (for GWASpoly and TASSEL)
# Uses three criteria: best GC, best replicability, and best significants
# PLINK also can produce info of more action models using options
#-----------------------------------------------------------
selectBestModel <- function (data, nBest, tool) {
	# Select main columns
	dataSNPs = data [,c("Marker","GC","MODEL","SCORE", "THRESHOLD", "DIFF")]; 

	# Order by nBest, DIFF, GC
	orderedDataSNPs = dataSNPs [order (dataSNPs$MODEL,-dataSNPs$DIFF),]; 

	for (N in c(200, 100, 50, nBest)) {
		# Reduce to groups of nBest
		groupsDataSNPs = Reduce (rbind, by(orderedDataSNPs, orderedDataSNPs["MODEL"], head, n=N)); 

		# Add Count of SNPs
		summ   = data.frame (add_count (groupsDataSNPs, Marker, sort=T)); 

		# Add count of significatives
		summSign   = data.frame (add_count (summ [summ$DIFF >0, ], "MODEL"))

		# Add count of shared SNPs
		summMdl = aggregate (x=summ$n, by=list(MODEL=summ$MODEL, GC=summ$GC), 	FUN=sum)

		# Summ differences
		#summMdlDiff = aggregate (x=summ$DIFF, by=list(MODEL=summ$MODEL, GC=summ$GC), 	FUN=sum)

		summMdlSign = cbind (summMdl, nSIGN=0)
		rownames (summMdlSign) = summMdlSign [,1]
		summMdlSign [as.character (summSign$MODEL),"nSIGN"] = ifelse(summSign$n==0, 0, summSign$n / sum (summSign$n))
			
		# Calculate best model score
		totalNs     = length (summMdlSign$MODEL) * N
		scoreGC     = 1 - abs (1-summMdlSign$GC)
		scoreShared = summMdlSign$x/totalNs  
		scoreSign   = summMdlSign$nSIGN 

		modelScore  = scoreGC + scoreShared + scoreSign
		summScores  = cbind (summMdl, scoreGC, scoreShared, scoreSign, score=modelScore)
		summScores  = summScores [order (summScores$score, summScores$MODEL, decreasing=T),]

		outFilename = paste0 ("report/model-scores-", tool, "-", sprintf ("%0.3d", N), ".csv")
		write.csv (summScores, file=outFilename, quote=F, row.names=F)
	}
	bestModel = summScores [1, "MODEL"]

	# Select SNPs for model and sort by DIFF
	dataModel = data [data[,"MODEL"] %in% bestModel,]
	dataModel = dataModel [order (-dataModel$DIFF),]

	return (dataModel)
}

#------------------------------------------------------------------------
# Create Venn diagram of common markers using info from summary table
#------------------------------------------------------------------------
markersVennDiagrams <- function (summaryTable, gwasModel, scoresType, outFile){
	# Params for figure shape and fonts
	WIDTH  = 7
	HEIGHT = 9
	CEXLABELS = 0.6
	CEXTITLES = 1.0

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
	v0 <- venn.diagram(x, height=2000, width=3000, alpha = 0.5, filename = NULL, # main=mainTitle,
						col = COLORS, cex=CEXLABELS, cat.cex=CEXTITLES, 
						margin=0.0, fill = COLORS)

	overlaps <- calculate.overlap(x)
	overlaps <- rev(overlaps)

	posOverlap = as.numeric (gsub ("a","", (names (overlaps))))
	for (i in 1:length(overlaps)){
		pos = posOverlap [i]
		v0[[pos+8]]$label <- paste(overlaps[[i]], collapse = "\n")
 	}


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
markersSummaryTable <- function (listOfResultsFile, gwasModel, nBest) {
	#files =  list.files(inputDir, pattern=sprintf ("(^(tool).*(%s).*[.](csv))", gwasModel), full.names=T)
	msgmsg ("Creating summary table...")
	summaryTable = data.frame ()

	tool=""
	for (resultsFile in listOfResultsFile) {
		data <- read.table (file=resultsFile, header=T)
		#if (nrow(data)>nBest) data=data [1:nBest,] 

		flagNewData = F
		if (grepl("GWASpoly", resultsFile)) {
			tool    = "GWASpoly"
			msgmsg (tool)
			data    = selectBestModel (data, nBest, tool)
			chrom   = data$Chrom
			pos	    = data$Position
			snps    = data$Marker
			flagNewData = T
		}else if (grepl ("PLINK", resultsFile)) {
			tool    = "PLINK"
			msgmsg (tool)
			chrom   = data$CHR
			pos	    = data$POS
			snps    = data$SNP
			flagNewData = T
		}else if (grepl ("TASSEL", resultsFile)) {
			tool    = "TASSEL"
			msgmsg (tool)
			data    = selectBestModel (data, nBest, tool)
			chrom   = data$Chr
			pos		= data$Pos
			snps    = data$Marker
			flagNewData = T
		}else if (grepl ("SHEsis", resultsFile)) {
			tool    = "SHEsis"
			msgmsg (tool)
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
hd1 = hd
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

	config = config::get (file=configFile) 

	configDF = data.frame (PARAMETER=character(), VALUE=character ())
	configDF = rbind  (configDF, data.frame (PARAMETER="Ploidy (4 or 2)", VALUE=toString (config$ploidy)))
	configDF = rbind  (configDF, data.frame (PARAMETER="Genotype filename", VALUE=toString (config$genotypeFile)))
	configDF = rbind  (configDF, data.frame (PARAMETER="Phenotype filename", VALUE=toString (config$phenotypeFile)))
	configDF = rbind  (configDF, data.frame (PARAMETER="Significance level (Genome-wide significance level)", VALUE=toString (config$significanceLevel)))
	configDF = rbind  (configDF, data.frame (PARAMETER="Correction method (Bonferroni or FDR)", VALUE=toString (config$correctionMethod)))
	configDF = rbind  (configDF, data.frame (PARAMETER="GWAS model (Full or Naive)", VALUE=toString (config$gwasModel)))
	configDF = rbind  (configDF, data.frame (PARAMETER="nBest (Number of best-ranked SNPs to be reported)", VALUE=toString (config$nBest)))
	configDF = rbind  (configDF, data.frame (PARAMETER="Filtering (TRUE or FALSE)", VALUE=toString (config$filtering)))
	configDF = rbind  (configDF, data.frame (PARAMETER="MIND Filter (Individual with missing genotype)", VALUE=toString (config$MIND)))
	configDF = rbind  (configDF, data.frame (PARAMETER="GENO Filter (SNPs with missing genotype)", VALUE=toString (config$GENO)))
	configDF = rbind  (configDF, data.frame (PARAMETER="MAF Filter (Minor allele frequency)", VALUE=toString (config$MAF)))
	configDF = rbind  (configDF, data.frame (PARAMETER="HWE Filter (Hardy-Weinberg test)", VALUE=toString (config$HWE)))
	configDF = rbind  (configDF, data.frame (PARAMETER="GWAS Tools", VALUE=toString (config$tools)))

	outName = paste0(outputDir, "/out-multiGWAS-inputParameters.tbl")
	write.table (file=outName, configDF, quote=F, sep="\t", row.names=F)
	config$workingDir = getwd ()
	return (config)
}

#-----------------------------------------------------------
# Create a chord diagram for SNPs vs Chromosomes from
# summary table of best scores
#-----------------------------------------------------------
createChordDiagramSharedSNPs <- function (scoresFile) {
	# >>>>> local function: matrix and colors 
	createChord <- function (matrixChord=NULL, colorsChord=NULL) {
		if (is.null (matrixChord)){
			plot.new()
			mtext ("No chord diagram (without shared SNPs)")
		}else {
			chordDiagram(matrixChord, annotationTrack = "grid", directional = -1, direction.type = c("arrows"), # c("diffHeight", "arrows"),
						 #link.arr.type = "big.arrow", 
						 grid.col = colorsChord,
						 preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(matrixChord))))))

			# we go back to the first track and customize sector labels
			circos.track(track.index = 1, panel.fun = function(x, y) {
					circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
					facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
			}, bg.border = NA) # here set bg.border to NA is important
		}
	} # >>>>> local function
	
	">>>>>> get shared SNPs <<<<<<>"
	getSharedSNPsFromFile <- function (scoresFile, N) {
		scores = read.table (file=scoresFile, header=T, sep="\t"); 
		summary = data.frame (add_count (scores, SNP, sort=T)); 
		sharedDups = summary [summary$n > 1,]
		shared = sharedDups [!duplicated (sharedDups$SNP),]
		return (shared)
	}
	">>>>>> get shared SNPs <<<<<<>"

	outFile = paste0(strsplit (scoresFile, split="[.]")[[1]][1], "-chordDiagram") 
	scores  = getSharedSNPsFromFile (scoresFile)

	# Check for shared SNPs
	matrixChord = NULL
	colorsChord = NULL
	if (nrow (scores) > 0) {
		# With shared SNPs
		tbl       = scores [,c("TOOL","CHROM","SNP")]
		tbl$CHROM = paste0 ("Chrom_", tbl$CHROM)

		# Group by TOOL and select 3 SNPs for each one
		#tblr = Reduce (rbind, by (tbl, tbl["TOOL"], head, n=2))
		tblm = tbl [,c(2,3)]
		chrs = sort (tblm [!duplicated (tblm [,1]), 1])
		snps = sort (as.character (tblm [!duplicated (tblm [,2]), 2]))

		# Create matrix Chroms X SNPs
		nChrs   = length (chrs)
		nSNPs   = length (snps)
		mat = as.data.frame (matrix (rep (0,nChrs*nSNPs),nrow=nChrs, ncol=nSNPs), stringAsFactor=F )
		rownames (mat) = chrs
		colnames (mat) = snps

		# Fill the matrix

		dmat = as.data.frame (mat)
		for (i in 1:nrow (tblm)) {
			chr = as.character (tblm [i, 1])
			snp = as.character (tblm [i, 2])
			dmat [chr,snp] = dmat [chr,snp] + 1
		}

		# Params for chordDiagram: mat and colors
		matrixChord = as.matrix (dmat) 

		colorsChrs = setNames (brewer.pal (n=nChrs+3, name="RdBu"), c(chrs, "chXX", "chYY", "chZZ"))
		colorsSNPs = setNames (rep ("grey", nSNPs), snps)
		colorsChord = c(colorsChrs, colorsSNPs)
	}

	# PDF
	pdf (file=paste0 (outFile, ".pdf"), width=7, height=7)
		createChord (matrixChord, colorsChord)
	dev.off()
	#PNG
	png (file=paste0 (outFile, ".png"), width=5, height=5	, units="in", res=90)
		createChord (matrixChord, colorsChord)
	dev.off()
}


#-------------------------------------------------------------
# Get alternate allele from a row of alleles
#-------------------------------------------------------------
#-------------------------------------------------------------
# Call to main function (first lines)
#-------------------------------------------------------------
#main ()

