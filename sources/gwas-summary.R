#!/usr/bin/Rscript

# INFO   : Create different summarized reports (tables, plots, Venn diagrams) for multiGWAS tool analysis
# AUTHOR : Luis Garreta (lgarreta@agrosavia.co)
# DATA   : feb/2020
# LOG: 
#	r5.2: Fixed unused factors in gwasResults
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
	source (paste0(HOME,"/sources/gwas-heatmap.R"))            # Module with functions to create summaries: tables and venn diagrams
	source (paste0(HOME,"/sources/gwas-preprocessing.R"))      # Module with functions to convert between different genotype formats 

	msgmsg ("Main...")

	inputDir      = "out/"
	genotypeFile  = "out/filtered-gwasp4-genotype.tbl"
	phenotypeFile = "out/filtered-gwasp4-phenotype.tbl"
	outputDir     = "report/"
	gwasModel     = "Full"
	nBest         = 10
	ploidy        = 4
	geneAction    = "all"

	createReports (inputDir, genotypeFile, phenotypeFile, 
				   ploidy, gwasModel, outputDir, nBest, geneAction)
}

#-------------------------------------------------------------
# Function to create different reports:
#	1- 1 table of best SNPs
#	2- 1 table of significative SNPs
#   3- 1 Venn diagram of best SNPs
#   4- 1 Venn diagram of significative SNPs
#	5- 1 multiplot of 4x4 manhattan and QQ plots
#-------------------------------------------------------------
#-------------------------------------------------------------
createReports <- function (inputDir, genotypeFile, phenotypeFile, ploidy, gwasModel, outputDir, nBest, geneAction) 
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
	if (length (listOfResultsFile)==0) {
		msgError ("WARNING: No result files for any tool. Check config file parameters (e.g. tools, geneAction, gwasModel)")
		quit ()
	}

	msgmsg ("Creating table with summary results...")
	snpTables = markersSummaryTable (listOfResultsFile, gwasModel, nBest, geneAction)

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
	op=markersManhattanPlots (listOfResultsFile, commonBest, snpTables, nBest, geneAction)
	dev.off()

	pdf (fileManhattanPlotPDF, width=11, height=15)
	op=markersManhattanPlots (listOfResultsFile, commonBest, snpTables, nBest, geneAction)
	par (op)
	dev.off()


	# Create heat maps
	msgmsg ("Creating SNP heatmaps for the best ranked SNPs...")
	genoNumericFilename = ACGTToNumericGenotypeFormat (genotypeFile, ploidy)

	createHeatmapForSNPList (outputDir, genotypeFile, genoNumericFilename, phenotypeFile, commonBest, ploidy)

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
markersManhattanPlots <- function (listOfResultsFile, commonBest, snpTables, nBest, geneAction) {
	#listOfResultsFile =  list.files(inputDir, pattern=sprintf ("(^(tool).*(%s).*[.](csv))", gwasModel), full.names=T)
	op <- par(mfrow = c(4,2), mar=c(3.5,3.5,3,1), oma=c(0,0,0,0), mgp = c(2.2,1,0))
	#op <- layout (matrix (c(1,1,2,3, 4:16), 4,4, byrow=T))#, widths=c(2,1), heights=c(2,1)))
	for (filename in listOfResultsFile) {
		data           = read.table (file=filename, header=T)
		data           = data [!is.na (data$P),]

		names          = unlist (strsplit (basename (filename), "[-|.]"))
		mainTitle      = paste0 (names[2],"-", names [3])

		if (grepl ("GWASpoly", filename)) {
			tool = "GWASpoly"
			data = selectBestModel (data, nBest, tool, geneAction)
			gwasResults = data.frame (SNP=data$Marker, CHR=data$Chrom, BP=data$Position, P=10^-data$SCORE, MODEL=data$MODEL)
		}
		else if (grepl ("PLINK", filename)) {
			tool = "PLINK"
			data = selectBestModel (data, nBest, tool, geneAction)
			gwasResults = data.frame (SNP=data$Marker, CHR=data$CHR, BP=data$POS, P=10^-data$SCORE, MODEL=data$MODEL)
		}
		else if (grepl ("TASSEL", filename)) {
			tool = "TASSEL"
			data = selectBestModel (data, nBest, tool, geneAction)
			gwasResults = data.frame (SNP=data$Marker, CHR=data$Chr, BP=data$Pos, P=10^-data$SCORE, MODEL=data$MODEL)
		}
		else if (grepl ("SHEsis", filename)) {
			tool = "SHEsis"
			gwasResults = data.frame (SNP=data$SNP, CHR=data$CHR, BP=data$POS, P=10^-data$SCORE, MODEL=data$MODEL)
		}

		# Remove old column factors
		gwasResults[] <- lapply(gwasResults, function(x) if(is.factor(x)) factor(x) else x)

		# Check for few rows
		nRows = nrow (data)
		nBest = ifelse (nBest>nRows, nRows, nBest )

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

		# Check if all chromosome names are numeric
		anyNonNumericChrom <- function (chrs) {
			suppressWarnings (any (is.na (as.numeric (chrs))))
		}

		# if non-numeric Chromosome names, convert to numeric using factors
		chrs = as.character (gwasResults$CHR)
		if  (anyNonNumericChrom (chrs)==TRUE) {
			msgmsg ("!!!Mapping chromosome names to numbers (see 'out-mapped-chromosome-names.csv') file...")
			chrs            = as.factor (chrs)
			levels (chrs)   = 1:length (levels (chrs))
			write.csv (data.frame (ORIGINAL_CHROM=gwasResults$CHR, NEW_CHROM=chrs), "out-mapped-chromosome-names.csv", quote=F, row.names=F)
			gwasResults$CHR = as.numeric (chrs)
		}

		msgmsg ("...Manhattan for", tool)
		manhattan(gwasResults,col = c("orange", "midnightblue"), highlight=sharedSNPs, annotatePval=bestThreshold, annotateTop=F,
				  suggestiveline=bestThresholdScore, genomewideline=signThresholdScore, main=mainTitle, logp=T, cex=2)

		text (x=0, y=signThresholdScore*0.92, "               Significants",, col="red", pos=4)
		text (x=0, y=bestThresholdScore*0.92, "Best",, col="blue", pos=4)

		#datax = calculateInflationFactor (-log10 (gwasResults$P))
		qqMGWAS (gwasResults, geneAction)
		#qq (gwasResults$P)
		#mtext (bquote(lambda[GC] == .(datax$delta)), side=3, line=-2, cex=0.7)
		#title (bquote(lambda[GC] == .(data$delta)))
	}
	#par (op)
	#dev.off()
	return (op)
}

#-------------------------------------------------------------
# QQ plot
#-------------------------------------------------------------
qqMGWAS <- function(gwasResults, geneAction) {
	models = levels (gwasResults$MODEL)
	#---- local fun -----
	qqValues <- function (pValues) {
		scores = -log10 (pValues)
		datax  = calculateInflationFactor (scores)
		n      = length(datax$scores)
		unif.p = -log10(ppoints(n))
		return (list (scores=scores, unif.p=unif.p, gc=datax$delta))
	}
	#---- local fun -----
	values = qqValues (gwasResults$P)
	xMax = max (values$unif.p)
	yMax = max (values$scores)

	#par(pty="s")
	#plot.new ()
	plot(1, type="n", xlim = c(0, xMax), ylim = c(0, yMax),
		 xlab=expression(paste("Expected -log"[10],"(p)",sep="")),
		 ylab=expression(paste("Observed -log"[10],"(p)",sep="")),
		 main=paste (models,collapse=" / "))
	lines(c(0,max(values$unif.p)),c(0,max(values$unif.p)),lty=2, col=c("red"))

	colors = c("black", "magenta3", "cyan")
	LEGEND = list()
	for (i in 1:length(models)) {
		pValues  = gwasResults [gwasResults$MODEL==models[i], c("P")]
		values=qqValues (pValues)
		points (values$unif.p, values$scores, pch=16, col=colors[i])
		GC = as.expression (bquote (.(models[i]) ~ lambda[GC] == .(values$gc)))
		LEGEND = append(LEGEND, GC) 
	}

	legend ("bottomright", lty=1, legend=LEGEND, col=colors, pch=16) 
}

#-----------------------------------------------------------
# Select best N SNPs from multiple action models (for GWASpoly and TASSEL)
# Uses three criteria: best GC, best replicability, and best significants
# PLINK also can produce info of more action models using options
#-----------------------------------------------------------
selectBestModel <- function (data, nBest, tool, geneAction) {
	if (geneAction!="all")
		return (data)

	# Select main columns
	dataSNPs = data [,c("Marker","GC","MODEL","SCORE", "THRESHOLD", "DIFF")]; 

	# Order by nBest, DIFF, GC
	orderedDataSNPs = dataSNPs [order (dataSNPs$MODEL,-dataSNPs$DIFF),]; 

	for (N in c(200, 100, 50, nBest)) {
		# Reduce to groups of nBest
		groupsDataSNPs = Reduce (rbind, by(orderedDataSNPs, orderedDataSNPs["MODEL"], head, n=N)); 

		# Add Count of SNPs between groups
		summ   = data.frame (add_count (groupsDataSNPs, Marker, sort=T, name="nSharedSNPs")); 

		# Add count of significatives
		summSign   = data.frame (add_count (summ [summ$DIFF >0, ], MODEL, name="nSign"))

		# Add count of shared SNPs
		summMdl = aggregate (x=summ$nSharedSNPs, by=list(MODEL=summ$MODEL, GC=summ$GC), 	FUN=sum)
		colnames (summMdl) = c("MODEL", "GC", "SHAREDNSPS")

		# Summ differences
		#summMdlDiff = aggregate (x=summ$DIFF, by=list(MODEL=summ$MODEL, GC=summ$GC), 	FUN=sum)

		# Add fraction of shared SNPs between all models
		summMdlSign = cbind (summMdl, nSIGN=0)
		rownames (summMdlSign) = summMdlSign [,1]
		summMdlSign [as.character (summSign$MODEL),"nSIGN"] = ifelse(summSign$nSharedSNPs==0, 0, summSign$nSharedSNPs / sum (summSign$nSharedSNPs))
			
		# Calculate best model score
		totalNs     = length (summMdlSign$MODEL) * N
		scoreGC     = 1 - abs (1-summMdlSign$GC)
		scoreShared = summMdlSign$SHAREDNSPS/totalNs  
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
	# Mark with "*" signficatives but if shared and no signifactive in one tool: not in the instersections
	#fs <- function (x) {if (x[10]=="TRUE") x[4]=paste0("*", x[4]);return(x)}
	#summaryTableMarks =as.data.frame (t(apply (summaryTable, 1, fs)))

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
markersSummaryTable <- function (listOfResultsFile, gwasModel, nBest, geneAction) {
	#files =  list.files(inputDir, pattern=sprintf ("(^(tool).*(%s).*[.](csv))", gwasModel), full.names=T)
	msgmsg ("Creating summary table...")
	summaryTable = data.frame (stringsAsFactors=F)

	tool=""
	for (resultsFile in listOfResultsFile) {
		data <- read.table (file=resultsFile, header=T, stringsAsFactors=F)
		#if (nrow(data)>nBest) data=data [1:nBest,] 

		flagNewData = F
		if (grepl("GWASpoly", resultsFile)) {
			tool    = "GWASpoly"
			msgmsg (tool)
			data    = selectBestModel (data, nBest, tool, geneAction)
			chrom   = data$Chrom
			pos	    = data$Position
			snps    = data$Marker
			flagNewData = T
		}else if (grepl ("PLINK", resultsFile)) {
			tool    = "PLINK"
			msgmsg (tool)
			data    = selectBestModel (data, nBest, tool, geneAction)
			chrom   = data$CHR
			pos	    = data$POS
			snps    = data$Marker
			flagNewData = T
		}else if (grepl ("TASSEL", resultsFile)) {
			tool    = "TASSEL"
			msgmsg (tool)
			data    = selectBestModel (data, nBest, tool, geneAction)
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
		pVal	= round (data$P, 6)
		pscores = round (data$SCORE, 4)
		tscores = round (data$THRESHOLD, 4)
		signf   = pscores >= tscores

		if (flagNewData==T) {
			dfm = data.frame (TOOL=tool, MODEL=model, GC=gcs, SNP=snps, CHROM=chrom, POSITION=pos, 
							  PVALUE=pVal, SCORE=pscores, THRESHOLD=tscores, SIGNIFICANCE=signf)
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

#-------------------------------------------------------------
# Get params from config file and define models according to ploidy
#-------------------------------------------------------------
## @knitr writeConfigurationParameters
writeConfigurationParameters <- function (inputDir, outputDir) {
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
						 grid.col = colorsChord, preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(matrixChord))))))

			# we go back to the first track and customize sector labels
			circos.track(track.index = 1, panel.fun = function(x, y) {
				xlim = get.cell.meta.data("xlim")
				xplot = get.cell.meta.data("xplot")
				ylim = get.cell.meta.data("ylim")
				sector.name = get.cell.meta.data("sector.index")							 

				#if(abs(xplot[2] - xplot[1]) < 20) {
				# Check if top or botton for Markers of Chromosomes
				if(abs(xplot[1]) < 180) 
					circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "red")
				 else 
					circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "blue")
				}, bg.border = NA) # here set bg.border to NA is important

				mtext ("Markers", side=3, col="red", cex=1.5)
				mtext ("Chromosomes", side=1, col="blue", cex=1.5)
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

		colorSNPs <- colorRampPalette(brewer.pal(8, "Set2"))(length(chrs))

		#colorsChrs = setNames (brewer.pal (n=nChrs+3, name="RdBu"), c(chrs, "chXX", "chYY", "chZZ"))
		colorsChrs = setNames (colorSNPs, c(chrs))
		colorsSNPs = setNames (rep ("grey", nSNPs), snps)
		colorsChord = c(colorsChrs, colorsSNPs)
	}

	funCreateChords <- function () {
		createChord (matrixChord, colorsChord)
	}

	# PDF
	pdf (file=paste0 (outFile, ".pdf"), width=7, height=7)
		funCreateChords ()
	dev.off()
	#PNG
	png (file=paste0 (outFile, ".png"), width=5, height=5	, units="in", res=90)
		funCreateChords ()
	dev.off()
}

#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
view <- function (data, m=10,n=10, tool="") {
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
msgmsg <- function (...) {
		messages = unlist (list (...))
		cat (">>>>", messages, "\n")
}

msgError <- function (...) {
		messages = unlist (list (...))
		cat (strrep("-", sum(sapply(messages, nchar))),"\n")
		cat (messages, "\n")
		cat (strrep("-", sum(sapply(messages, nchar))),"\n")
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
# Call to main function (first lines)
#-------------------------------------------------------------

#source ("lglib06.R")
#main ()

