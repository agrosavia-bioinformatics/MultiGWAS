
#-------------------------------------------------------------
# SHEsis tool
# Shesis uses a low significance level of 0.01 as their pvalues are inflates
#-------------------------------------------------------------
runToolShesis <- function (params) 
{
	geneAction = params$geneAction
	if (geneAction != "all" & geneAction != "additive")
		return ()

	geneAction = "additive"
	msgmsg ("Running SHEsis GWAS...")

	outFile      = paste0 ("out/tool-SHEsis-scores-", params$gwasModel)
	scoresFile   = paste0 (outFile,".csv")

	if (params$gwasModel == "Naive") {
		gwaspToShesisGenoPheno (params$genotypeFile, params$phenotypeFile, params$ploidy, params$traitType)
	} else if (params$gwasModel == "Full") {
		# Apply kinship and filter individuals
		inGeno  = "out/filtered-plink-genotype"       # Uses plink file
		kinFile = paste0 (inGeno,"-kinship-shesis")
		cmm     = sprintf ("%s/sources/scripts/script-kinship-plink2.sh %s %s", HOME, inGeno, kinFile)
		runCommand (cmm, "log-kinship.log")

		# Filter geno/pheno to individuals, write out, and convert to SHEsis format
		kinshipIndividuals  = read.table (file=paste0(kinFile, ".ped"), sep="\t", check.names=F)
		individuals         = as.character (kinshipIndividuals [,1])

		genotype             = read.csv (file=params$genotypeFile, header=T, check.names=F)
		phenotype            = read.csv (file=params$phenotypeFile, header=T, check.names=F)
		genotypeKinship      = genotype [, c(colnames (genotype)[1:3], individuals)]
		phenotypeKinship     = phenotype [phenotype [,1] %in% individuals,]
		genotypeFileKinship  = addLabel (params$genotypeFile,  "kinship-shesis")
		phenotypeFileKinship = addLabel (params$phenotypeFile, "kinship-shesis")
		write.table (file=genotypeFileKinship,  genotypeKinship,  row.names=F, quote=F, sep=",")
		write.table (file=phenotypeFileKinship, phenotypeKinship, row.names=F, quote=F, sep=",")

		gwaspToShesisGenoPheno (genotypeFileKinship, phenotypeFileKinship, params$ploidy, params$traitType)
	}

	runShesisCommand (params$traitType, outFile, params, geneAction)
	msg ("... Ending SHEsis")
}

#-------------------------------------------------------------
#-------------------------------------------------------------
runShesisCommand <- function (traitType, outFile, params, geneAction) {
	inGenoPheno  = "out/filtered-shesis-genopheno.tbl"
	inMarkers    = "out/filtered-shesis-markernames.tbl"
	flagQTL      = ifelse (traitType=="quantitative", "--qtl", "")
	scoresFile   = paste0 (outFile, ".csv")

	cmm=sprintf ("%s/sources/scripts/script-shesis-associations-qtl.sh %s %s %s %s %s", HOME, inGenoPheno, params$ploidy, inMarkers, outFile, flagQTL)
	runCommand (cmm, "log-SHEsis.log")	

	if (traitType=="quantitative")
		resultsAll = createTableFromQuantitativeResults (outFile, params, geneAction)
	else # case-control
		resultsAll = createTableFromBinaryResults (outFile, params, geneAction)

	write.table (file=scoresFile, resultsAll, row.names=F, quote=F, sep="\t")
	return (resultsAll)
}

#-------------------------------------------------------------
# Create table from binart .txt results file
#-------------------------------------------------------------
createTableFromBinaryResults <- function (outFile, params, geneAction) {
	resultsFile = paste0 (outFile, ".txt")

	conn = file (resultsFile, open="r")
	lines = readLines (conn)
	close (conn)
	n  = length (lines) 

	inAllele=F; genotypeList = list (); allelesList  = list (); currentMarker = NULL
	for (i in 1:n)  {
		li = lines [i]
		if (grepl ("Allele", li)) {
			inAllele=T
			marker        = strsplit (li, "(", fixed=T)[[1]][1]
			currentMarker = c(SNP=marker)
		}else if (grepl ("Genotype", li)) {
			inAllele=F
			marker        = strsplit (li, "(", fixed=T)[[1]][1]
			currentMarker = c(SNP=marker)
		}else if (grepl ( "Pearson", li)) {
			p = strsplit (li, " ")
			currentMarker = c (currentMarker, pPearson=p[[1]][length(p[[1]])])
			if (inAllele)
				allelesList = append (allelesList, currentMarker)
			else
				genotypeList = append (genotypeList, currentMarker)
		}
	}

	columns = c("SNP", "pPearson")
	allelesMat   = matrix (allelesList, ncol=2, byrow=T)
	colnames (allelesMat) = columns
	write.csv (allelesMat, "shesis-alleles-pvalues.csv", quote=F, row.names=F)

	results   = matrix (genotypeList, ncol=2, byrow=T)
	colnames (results) = columns
	write.csv (results, "shesis-genotype-pvalues.csv", quote=F, row.names=F)
	results = read.csv ("shesis-genotype-pvalues.csv")

# LG: Added 1e-10 to avoid "inf" values in scores
	#pValues  = results[,"P.value"] + 1e-10

	pValues  = results[,"pPearson"] 
	print (pValues)
	m        = length (pValues)
	adj       = adjustPValues (0.01, pValues, params$correctionMethod)
	message ("Formating 0 ...")
	pValues   = adj$pValues
	message ("Formating 1 ...")
	threshold = adj$threshold
	message ("Formating 2 ...")
	scores    = -log10 (pValues)

	# Set Columns
	map <- read.table (file="out/map.tbl", sep="\t", check.names=F)
	rownames (map) = map [,1]
	SNP <- as.character (results$SNP)
	CHR <- map [SNP, 2]
	POS <- map [SNP, 3]
	GC  = calculateInflationFactor (scores)

	model = geneAction
	resultsAll <- data.frame (MODEL=model, GC=GC$delta, SNP, CHR, POS, P=pValues, SCORE=round (scores,6), THRESHOLD=round (threshold,6), DIFF=round (scores-threshold, 6), results)
	resultsAll <- resultsAll [order (resultsAll$DIFF, decreasing=T),]

	return (resultsAll)
}

#-------------------------------------------------------------
# Create table from quantitative .txt file 
#-------------------------------------------------------------
createTableFromQuantitativeResults <- function (outFile, params, geneAction) {
	resultsFile = paste0 (outFile, ".txt")

	results  = read.table (file=resultsFile, header=T, sep="\t", check.names=T) # TRUE as SHEsis colnames have spaces

	# LG: Added 1e-10 to avoid "inf" values in scores
	#pValues  = results[,"P.value"] + 1e-10
	pValues  = results[,"P.value"] 
	m        = length (pValues)
	adj       = adjustPValues (0.01, pValues, params$correctionMethod)
	pValues   = adj$pValues
	threshold = adj$threshold
	scores    = -log10 (pValues)

	# Set Columns
	map <- read.table (file="out/map.tbl", sep="\t", check.names=F)
	rownames (map) = map [,1]
	SNP <- as.character (results$SNP)
	CHR <- map [SNP, 2]
	POS <- map [SNP, 3]
	GC  = calculateInflationFactor (scores)

	model = geneAction
	resultsAll <- data.frame (MODEL=model, GC=GC$delta, SNP, CHR, POS, P=pValues, SCORE=round (scores,6), THRESHOLD=round (threshold,6), DIFF=round (scores-threshold, 6), results)
	resultsAll <- resultsAll [order (resultsAll$DIFF, decreasing=T),]

	return (resultsAll)
}

#----------------------------------------------------------
# Transform table genotype to SHEsis genotype format
#----------------------------------------------------------
gwaspToShesisGenoPheno <- function (genotypeFile, phenotypeFile, ploidy, traitType) 
{
	msgmsg ("    >>>> Writting gwasp to shesis genopheno...")
	sep <- function (allele) {
		s="";
		for (i in 1:ploidy) s=paste0(s, substr (allele,start=i,stop=i)," ");
		#s = paste (strsplit (allele, "")[[1]], collapse=" ")
		return (s)
	}
	geno    = read.csv (file=genotypeFile, stringsAsFactors=F, check.names=F)
	pheno   = read.csv (file=phenotypeFile, stringsAsFactors=F, check.names=F)
	if (traitType=="case-control") {
		pheno [pheno[,2]==1,2]=2
		pheno [pheno[,2]==0,2]=1
	}

	rownames (pheno) <- pheno [,1]
	map        <- geno  [,c(1,2,3)]    # Get SNP, Cromosome, Position
	rownames (geno)  <- map [,1] 

	alleles    <- geno[,-c(1,2,3)]
	alleles [is.na(alleles)] = paste (rep ("0", ploidy), collapse="") # NAs as "00" or "0000"

	allelesMat <- t(sapply (alleles, sep))

	samples         = rownames (allelesMat)
	pheno           = pheno [samples,]
	pheno [,2]      = impute.mode (pheno [,2])
	genoPhenoShesis = data.frame (Sample=pheno[,1], Trait=pheno[,2],  allelesMat)

	msgmsg ("    >>>> Writing shesis genopheno...")
	outFile = "out/filtered-shesis-genopheno.tbl"
	write.table (file=outFile, genoPhenoShesis, quote=F,row.names=F,col.names=F, sep="\t")

	msgmsg ("    >>>> Writing shesis marker names...")
	outFile = "out/filtered-shesis-markernames.tbl"
	write.table (file=outFile, map[,1], quote=F,row.names=F,col.names=F, sep="\t")
}

