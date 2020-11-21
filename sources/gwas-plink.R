#-------------------------------------------------------------
# Plink tool
# Plink also has gene action models: ADD, DOM, ADD+DOM (GNR)
# https://zzz.bwh.harvard.edu/plink/anal.shtml:
# The TEST column is by default ADD meaning the additive effects of allele dosage. Adding the option
#     --genotypic
# will generate file which will have two extra tests per SNP, corresponding to two extra rows: 
# DOMDEV and GENO_2DF which represent a separate test of the dominance component or a 2 df joint 
# test of both additive and dominance (i.e. corresponding the the general, genotypic model in the
# --model command). Unlike the dominance model is the --model, DOMDEV refers to a variable coded 
# 0,1,0 for the three genotypes AA,Aa,aa, i.e. representing the dominance deviation from additivity, 
# rather specifying that a particular allele is dominant or recessive. That is, the DOMDEV term is 
# fitted jointly with the ADD term in a single model.
#-------------------------------------------------------------

runToolPlink <- function (params) {
	msgmsg ("Running Plink GWAS...")

	outFile = paste0 ("out/tool-PLINK-scores-", params$geneAction,"-", params$gwasModel)

	if (params$geneAction=="additive")
		scores = runPlinkCommand (params, '', outFile)
	else if (params$geneAction=="general")
		scores = runPlinkCommand (params, "genotypic", outFile)
	else if (params$geneAction=="dominant") {
		recScores = runPlinkCommand (params, "recessive", outFile)
		domScores = runPlinkCommand (params, "dominant", outFile)
		scores     = rbind (recScores, domScores)
	}else if (params$geneAction=="all") {
		addScores = runPlinkCommand (params, '', outFile)
		gnrScores = runPlinkCommand (params, "genotypic", outFile)
		recScores = runPlinkCommand (params, "recessive", outFile)
		domScores = runPlinkCommand (params, "dominant", outFile)
		scores    = rbind (addScores, gnrScores, recScores, domScores)
	}
	scores     = scores [order (scores$DIFF, decreasing=T),]
	scoresFile = paste0 (outFile, ".csv")
	colnames (scores) = gsub ("SNP", "Marker", colnames (scores))
	write.table (scores, scoresFile, row.names=F, quote=F, sep="\t")
	msg ("... Ending PLINK")
}

#-------------------------------------------------------------
#-------------------------------------------------------------
runPlinkCommand <- function (params, geneAction, outFile) 
{
	model   = params$gwasModel
	inGeno           = "out/filtered-plink-genotype"       # Only prefix for Plink
	inPheno          = "out/filtered-plink-phenotype.tbl"  
	outPlinkLinear   = paste0 (outFile,".TRAIT.assoc")

	# Type of trait: quantitative (e.g 1.0023 or case-control (e.g. 0/1)
	if (params$traitType=="quantitative") {
		outPlinkLinear   = paste0 (outPlinkLinear, ".linear")
		flagTrait = ""
	}else {
		outPlinkLinear   = paste0 (outPlinkLinear, ".logistic")
		flagTrait = "--1"
	}

	# Naive or FUll GWAS model
	if (model=="Naive") { 
		cmm=sprintf ("%s/sources/scripts/script-plink-NaiveModel.sh %s %s %s", HOME, inGeno, inPheno, outFile)
		runCommand (cmm, "log-Plink.log")

	}else if (model=="Full") {
		# FIRST: kinship filtering of .ped file
		outKinFile = paste0 (inGeno, "-", geneAction, "-kinship-plink")
		cmm=sprintf ("%s/sources/scripts/script-kinship-plink2.sh %s %s", HOME, inGeno, outKinFile)
		runCommand (cmm, "log-kinship.log")

		# SECOND: Structure by PCs 
		cmm=sprintf ("%s/sources/scripts/script-plink-FullModel.sh %s '%s' '%s' %s %s", HOME, outKinFile, flagTrait, inPheno, outFile, geneAction)
		runCommand (cmm, "log-Plink.log")
	}

	# Read from the secondary output file from PLINK
	resultsLinearAll = read.table (file=outPlinkLinear,  header=T, check.names=F) 
	resultsLinear    = resultsLinearAll [!duplicated (resultsLinearAll$SNP),]
	resultsLinear    = resultsLinear [!is.na (resultsLinear$P),]

	pValues     = resultsLinear$P
	chromosomes = resultsLinear$CHR
	positions   = resultsLinear$BP

	# Adjust scores, threshold using Bonferroni or FDR
	adj       = adjustPValues (params$significanceLevel, pValues, params$correctionMethod)
	pValues   = adj$pValues
	threshold = round (adj$threshold, 6)
	scores    = round (-log10 (pValues), 6)
	gcs       = calculateInflationFactor (scores)
	diffs     = round (scores-threshold, 6)

	model      = ifelse (geneAction=='', "additive", geneAction)  # But Plink offers other tow gene action models: "ADD", "DOM", "GNR" (ADD+COM)
	resultsAll = cbind (MODEL=model, GC=gcs$delta, CHR=chromosomes, POS=positions, P=pValues, 
						SCORE=scores, THRESHOLD=threshold, DIFF=diffs, resultsLinear)
	return (resultsAll)
}
