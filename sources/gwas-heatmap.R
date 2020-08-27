#!/usr/bin/Rscript
# LOG: 
#     Jul 11: Removed NAs with complete.cases
#     Jun 06: Fixed, not include NAs in both geno and pheno

suppressMessages (library (gplots)) 
#----------------------------------------------------------
# Main
#----------------------------------------------------------
main <- function () {
	args = commandArgs (trailingOnly=T)

	genoFileACGT <- "geno-ACGT.csv"
	genoFileNUM  <- "geno-NUM.csv"
	phenoFile    <- "pheno.csv"
	snpList      <- c ("c1_8019", "c2_21314")         # SNPs a visualizar
	createHeatmapForSNPList ("./",genoFileACGT, genoFileNUM, phenoFile, snpList) 
}

#----------------------------------------------------------
# create SNP profiles for a list of SNPS
#----------------------------------------------------------
createHeatmapForSNPList <- function (outputDir, genoFileACGT, genoFileNUM, phenoFile, snpList, ploidy) {
	#msg ("Genotype ACGT : ", genoFileACGT)
	#msg ("Genotype NUM  : ", genoFileNUM)
	#msg ("Phenotype     : ", phenoFile)
	#message ("SNPs:         : ", snpList)

	outName = paste0 (outputDir, "/out-SNPProfile") 

	pdfHeatMap <- function (snp) {
		message  ("    >>>> Heatmap for snp: ", snp)
		createHeatmapForSNP (outputDir, genoFileACGT, genoFileNUM, phenoFile, snp, ploidy)
	}

	NCORES = detectCores ()
	for (s in snpList)
		pdfHeatMap (s)
	#res=mclapply (snpList, pdfHeatMap, mc.cores=NCORES)
}

#----------------------------------------------------------
# create a SNP profile for a snpId
#----------------------------------------------------------
createHeatmapForSNP <- function (outputDir, genoFileACGT, genoFileNUM, phenoFile, snpId, ploidy) {
	genotypeACGT    <- read.csv (genoFileACGT, na.strings = "NA", dec = ".", strip.white = TRUE, check.names=F )
	genotypeNUMERIC <- read.csv (genoFileNUM,  na.strings = "NA", dec = ".", strip.white = TRUE, check.names=F)
	phenotype       <- read.csv (phenoFile,    na.strings = "NA", dec = ".", strip.white = TRUE, check.names=F)

	# Get names and values for phenotype
	phenoNames  = phenotype [,1]
	phenoValues = phenotype [,2]
	trait       = colnames (phenotype) [2]


	# Get samples for marker in both matrices: Numeric and ACGT
	samplesMarker          = t(genotypeNUMERIC[genotypeNUMERIC[,1]==snpId,])
	samplesMarkerMatrixNUM = as.matrix(samplesMarker[-1:-3,])

	samplesMarker           = t(genotypeACGT[genotypeACGT[,1]==snpId,])
	samplesMarkerMatrixACGT = as.matrix(samplesMarker[-1:-3,])

	# 
	one.hot<-function(sqnce, alphabet){
	  y<-unlist(strsplit(sqnce,""))
	  sapply(y,function(x){match(alphabet,x,nomatch=0)})
	}

	if (ploidy==4) alphabet_AA = c("0","1","2","3","4")
	else           alphabet_AA = c("0","1","2") 

	c           <- t(one.hot (samplesMarkerMatrixNUM, alphabet = alphabet_AA))
	names(c)    <- names(samplesMarkerMatrixNUM)
	e           <- as.data.frame(c,row.names = samplesMarkerMatrixNUM[1,])
	e$Name      <- row.names(samplesMarkerMatrixNUM)
	genoxfeno_2 <- subset(e,e$Name%in%phenoNames)

	if (ploidy==4)
		genoxfeno_3 <- genoxfeno_2[,1:5]*phenoValues
	else
		genoxfeno_3 <- genoxfeno_2[,1:3]*phenoValues

	# Remove rows with NAs from genoxfeno_3 
	g3NAs                 = genoxfeno_3
	rownames (g3NAs)      = phenoNames 
	g3noNAs               = g3NAs[complete.cases(g3NAs),]
	genoxfenov4           = as.matrix(g3noNAs)
	rownames(genoxfenov4) = rownames (g3noNAs)

	data_geno_NUM_ACGT<-data.frame(samplesMarkerMatrixNUM,samplesMarkerMatrixACGT)
	
	gen_0 <- subset (data_geno_NUM_ACGT, samplesMarkerMatrixNUM=="0",select = samplesMarkerMatrixACGT)
	gen_1 <- subset (data_geno_NUM_ACGT, samplesMarkerMatrixNUM=="1",select = samplesMarkerMatrixACGT)
	gen_2 <- subset (data_geno_NUM_ACGT, samplesMarkerMatrixNUM=="2",select = samplesMarkerMatrixACGT)

	g0    <- as.character (gen_0$samplesMarkerMatrixACGT[1])    
	g1    <- as.character (gen_1$samplesMarkerMatrixACGT[1])
	g2    <- as.character (gen_2$samplesMarkerMatrixACGT[1])

	if (ploidy==4) {
		gen_3 <- subset (data_geno_NUM_ACGT, samplesMarkerMatrixNUM=="3",select = samplesMarkerMatrixACGT)
		gen_4 <- subset (data_geno_NUM_ACGT, samplesMarkerMatrixNUM=="4",select = samplesMarkerMatrixACGT)

		g3    <- as.character (gen_3$samplesMarkerMatrixACGT[1])
		g4    <- as.character (gen_4$samplesMarkerMatrixACGT[1])
		colnames(genoxfenov4)<-c(g0,g1,g2,g3,g4)
	}else
		colnames(genoxfenov4)<-c(g0,g1,g2)


	my_palette <- colorRampPalette(c("white", "blue", "darkblue", "darkred"))(n = 100)
	lmat <- rbind(c(5,4), c(2,3), c(2,1))
	lhei <- c(12,0.1,32)
	lwid <- c(3,9)

	myplot<- function() {
	  oldpar <- par("mar")
	  hist(phenoValues, main = "Histogram", xlab=trait)
	}

	# Function that call heatmap.2 (only to save code)
	fun_heatmap <- function () {
		hmap = heatmap.2(genoxfenov4,dendrogram = "row",reorderfun=function(snpId, w) reorder(snpId, w, agglo.FUN = mean),
				  Colv=FALSE,adjCol = c(NA,0),key=TRUE,srtCol=360,
				  col=my_palette,cexCol = 1,lmat=lmat,lhei=lhei, lwid=lwid, extrafun=myplot, 
				  key.xlab=paste0 ("Value of ", trait), xlab=snpId, key.title="Color Key")
	}

	outName = paste0 (outputDir, "/out-SNPProfile") 
	pdf (paste0(outName,"-", snpId, ".pdf"), width=7, height=7)
		fun_heatmap ()
	dev.off()

	png (paste0(outName,"-",snpId, ".png"), width=7, height=7, units="in", res=72)
		fun_heatmap ()
	dev.off()
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
		messages = unlist (list (...))
		cat (">>>>", messages, "\n")
}

#----------------------------------------------------------
#----------------------------------------------------------
#main()
