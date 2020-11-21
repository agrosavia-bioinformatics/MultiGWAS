#!/usr/bin/Rscript
DEBUG = F
options (width=300)
if (DEBUG) options (warn=2)
#source ("lglib06.R")

# INFO  : Tool for running GWAS integratind four GWAS tools: 
#         GWASpoly and SHEsis for polyploids species, and Plink and Tassel for diploids.
# AUTHOR: Luis Garreta (lgarreta@agrosavia.co) 
# DATE  : 12/feb/2020
# LOGS  :   
	# r1.5: Using VCF files
	# r1.3: Added column gene action model to tables of results by tool
	# r1.2: Separate SHEsis/Plink Kinship. Reduced code (runGWAStools, ped, bed)

#-------------------------------------------------------------
# Return string with usage instructions
#-------------------------------------------------------------
usageInstructions <- function () {
	USAGE="USAGE: multiGWAS <config file>"
	return (USAGE)
}

#-------------------------------------------------------------
# Main for multi traits
#-------------------------------------------------------------
main <- function () {
	message ("MultiGWAS 1.0")
	args = commandArgs(trailingOnly = TRUE)

	if (length (args) < 1) 
		stop (usageInstructions())
	else if (substring (args [1],1,2)=="--")
		processCommandArguments (args)
	else
		processConfigFile (args)
}

#-------------------------------------------------------------
# Define global variables, load packages, and load sources
#-------------------------------------------------------------
initGlobalEnvironment <- function () {
	msg ("Loading libraries and setting globals...")
	# Get enviornment GWAS variable 
	HOME <<- Sys.getenv ("MULTIGWAS_HOME")
	.libPaths (paste0(HOME, "/opt/Rlibs"))

	# Load packages
	suppressMessages (library (GWASpoly)) #
	suppressMessages (library (parallel)) #
	suppressMessages (library (config))  # For read config file
	NCORES <<- ifelse (DEBUG==T, 1, detectCores ())

	# New class for gwaspoly
	setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")

	# Load sources
	source (paste0 (HOME, "/sources/gwas-preprocessing.R"))      # Module with functions to convert between different genotype formats 
	source (paste0 (HOME, "/sources/gwas-summary.R"))            # Module with functions to create summaries: tables and venn diagrams
	source (paste0 (HOME, "/sources/gwas-heatmap.R"))            # Module with functions to create heatmaps for shared SNPs
	source (paste0 (HOME, "/sources/gwas-gwaspoly.R"))           # Module with gwaspoly functions
	source (paste0 (HOME, "/sources/gwas-plink.R"))              # Module with plink functions
	source (paste0 (HOME, "/sources/gwas-tassel.R"))             # Module with tassel functions
	source (paste0 (HOME, "/sources/gwas-shesis.R"))             # Module with shesis functions
	source (paste0 (HOME, "/sources/gwas-lib.R"))             # Module with shesis functions
}

#-------------------------------------------------------------
# Process command line arguments
#-------------------------------------------------------------
processCommandArguments <- function (args) {
	if (args [1]=="--fitpoly")
		commandArgsFitpoly (args)
}

#-------------------------------------------------------------
# Process config file and run multiGWAS
#-------------------------------------------------------------
processConfigFile <- function (args) {
	msg ("Processing config file...")
	initGlobalEnvironment ()

	configFile = args [1]
	if (file.exists (configFile)==F) 
		stop ("Configuration file not found")

	# Read and check config file arguments
	msg ("Reading configuration file...")
	config     = getMainConfig (configFile)

	mainDir = getwd ()
	for (traitConfigName in config$traitConfigList) {
		setwd (mainDir)
		mainSingleTrait (traitConfigName)
	}
}

#-------------------------------------------------------------
# Create individual phenotype files from a multitrait phenotype file
#-------------------------------------------------------------
createTraitConfigFiles <- function (phenotype, traitName, configFile, config) {
	# Create phenotype file for trait
	phenotypeTrait     = phenotype [, c(colnames(phenotype)[1], traitName)]
	phenotypeTraitPath = paste0 (traitName,".csv")
	write.csv (phenotypeTrait, phenotypeTraitPath, quote=F, row.names=F)

	# Create config file for trait
	config$phenotypeFile = basename (phenotypeTraitPath)
	configLines = c("default:")
	for (n in names (config))
		configLines = c (configLines, paste0 ("  ", n, " : ", config[n]))

	configLines      = c (configLines, paste0 ("  outputDir : ", traitName))
	traitConfigName  = paste0 (traitName, ".config")
	writeLines (configLines, traitConfigName)

	return (traitConfigName)
}

#-------------------------------------------------------------
# Main for a single trait
#-------------------------------------------------------------
mainSingleTrait <- function (traitConfigName) {
	# Copy files to working dirs and create output dirs
	config     = getTraitConfig (traitConfigName)

	# Read, filter, and check phenotype and genotype
	msg ("Preprocessing genomic data (Filtering and Formating data)...")

	data <- genoPhenoMapProcessing (config$genotypeFile, config$genotypeType,
									config$phenotypeFile, config$mapFile, config)

	config$genotypeFile  = data$genotypeFile
	config$phenotypeFile = data$phenotypeFile
	config$trait         = data$trait

	# Run the four tools in parallel
	runGWASTools (config)

	# Create reports
	msg ("Creating reports (Table, Venn diagrams, Manhattan&QQ plots, SNP profiles)...")
	createReports (config$outputDir, config$genotypeFile, config$phenotypeFile, config$ploidy,
				   config$gwasModel, config$reportDir, nBest=as.integer (config$nBest), config$geneAction)

	# Move out files to output dir
	msg ("Moving files to output folders...")
	moveOutFiles (config$outputDir, config$reportDir)
}
#-------------------------------------------------------------
# Read configuration parameters and set working directories
#   - Create dir, if it exists, it is renamed as old-XXXX
#   - Copy and make links of geno/pheno to output dirs
#-------------------------------------------------------------
getTraitConfig <- function (traitConfigName) {
	msg ("Processing config file: ", traitConfigName)
	config     = config::get (file=traitConfigName, config="advanced") 
	traitDir   = config$outputDir
	outDir     = paste0 (traitDir, "/out")


	# Copy files to trait dir and out dir
	outDirs    = c(traitDir, outDir)
	for (dir in outDirs) {
		createDir (dir)
		runCommand (sprintf ("cp %s %s", traitConfigName, dir))
		runCommand (sprintf ("cp %s %s", config$genotypeFile, dir ))
		runCommand (sprintf ("cp %s %s", config$phenotypeFile, dir))
		if (config$genotypeType %in% c("kmatrix", "fitpoly"))
			runCommand (sprintf ("cp %s %s", config$mapFile, dir))
	}
	# Change to the working dir and set dirs in config
	setwd (traitDir)
	config$outputDir <- "out/"
	config$reportDir <- "report/"

	return (config)
}

#-------------------------------------------------------------
# Get params from config file and define models according to ploidy
#-------------------------------------------------------------
getMainConfig <- function (configFile) {
	config     = config::get (file=configFile, config="advanced") 
	config$configFilename = configFile

	# Still unimplemented "all" option in some tools for (e.g PLINK)
	if (is.null (config$geneAction))
		config$geneAction = "additive"
	if (is.null (config$traitType))
		config$traitType = "quantitative"

	# Check input files
	if (!file.exists (config$genotypeFile) | !file.exists (config$phenotypeFile)) {
		errorMessage = "Input files (genotype or phenotype) not found!!"
		errorMessage = paste0 (errorMessage, "\n\t Genotype: ", config$genotypeFile)
		errorMessage = paste0 (errorMessage, "\n\t Phenotype: ", config$phenotype)
		stop (errorMessage, call.=T)
	}
	if (config$genotypeType=="kmatrix" && !file.exists (config$mapFile)) {
		errorMessage = "map file not found!!"
		errorMessage = paste0 (errorMessage, "\n\t Map file: ", config$mapFile)
		stop (errorMessage)
	}

	# Print config file
	msgmsg ("------------------------------------------------")
	msgmsg ("Summary of configuration parameters:")
	msgmsg ("------------------------------------------------")
	msgmsg ("Ploidy                 : ", config$ploidy) 
	msgmsg ("Genotype filename      : ", config$genotypeFile) 
	msgmsg ("Phenotype filename     : ", config$phenotypeFile) 
	msgmsg ("Significance level     : ", config$significanceLevel) 
	msgmsg ("Correction method      : ", config$correctionMethod) 
	msgmsg ("GwAS model             : ", config$gwasModel) 
	msgmsg ("Number of repored SNPs : ", config$nBest) 
	msgmsg ("Filtering              : ", config$filtering) 
	msgmsg ("MIND                   : ", config$MIND) 
	msgmsg ("GENO                   : ", config$GENO) 
	msgmsg ("MAF                    : ", config$MAF) 
	msgmsg ("HWE                    : ", config$HWE) 
	msgmsg ("Tools                  : ", config$tools) 
	msgmsg ("------------------------------------------------")
	msgmsg ("Gene action model      : ", config$geneAction) 
	msgmsg ("Trait type             : ", config$traitType) 
	msgmsg ("------------------------------------------------")

	# Create output dir for this project
	outDir   = paste0 ("out-", strsplit (configFile, split="[.]") [[1]][1])
	createDir (outDir)
	runCommand(sprintf ("cp %s %s", config$genotypeFile, outDir))
	runCommand(sprintf ("cp %s %s", config$phenotypeFile, outDir))
	runCommand(sprintf ("cp %s %s", config$mapFile, outDir))

	# Change to the working dir of the main projet 
	setwd (outDir)
	
	# Create config files for each trait
	phenotype = read.csv (config$phenotypeFile, check.names=F)
	traitList = colnames (phenotype)[-1]
	traitConfigList = c()
	for (traitName in traitList) {
		traitConfig     = createTraitConfigFiles (phenotype, traitName, configFile, config)
		traitConfigList = c (traitConfigList, traitConfig)
	}

	config$traitConfigList = traitConfigList 
	return (config)
}

#-------------------------------------------------------------
# Used to run in parallel the other functions
#-------------------------------------------------------------
runGWASTools <- function (config) {
	runOneTool <- function (tool, config) {
		if      (tool=="gwaspoly") runToolGwaspoly (config)
		else if (tool=="plink")    runToolPlink (config)
		else if (tool=="shesis")   runToolShesis (config)
		else if (tool=="tassel")   runToolTassel (config)
		else                       stop ("Tool not supported")
	}

	# A string containing the names of the tools to run (e.g. "GWASpoly SHEsis PLINK TASSEL")
	#config$tools = "Plink Shesis"
	config$tools = strsplit(tolower (config$tools) ,split=" ")[[1]]
	msg ("Preparing to execute in parallel the GWAS tools:") 
	for (i in 1:length(config$tools)) 
		msgmsg ("Running ", config$tools [i])

	mclapply (config$tools, runOneTool, config, mc.cores=NCORES, mc.silent=T)
}


#-------------------------------------------------------------
# Move output files to specific directories
#-------------------------------------------------------------
moveOutFiles <- function (outputDir, reportDir) 
{
	system (sprintf ("cp %s/tool*csv %s &> /dev/null", outputDir, reportDir))
	system (sprintf ("mv %s/out*pdf %s > /dev/null 2>&1", outputDir, reportDir))
	system ("mkdir logs")
	system ("mv *.log* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv *.errors* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv *PCs* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv ../*log* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
	system ("mv ../*errors* logs > /dev/null 2>&1", ignore.stdout=T, ignore.stderr=T)
}

#-------------------------------------------------------------
# Filters the genotype by different quality control filters
# Read and check files, sample samples
# Convert geno an pheno to other tool formats 
#-------------------------------------------------------------
genoPhenoMapProcessing <- function (genotypeFile, genotypeType, phenotypeFile, mapFile, config) {
	# Check for VCF, GWASpoly, k-matrix, or fitPoly.
	# and convert to GWASpoly format
	genotypeFile = genotypeProcessing (genotypeFile, genotypeType, mapFile, config$ploidy)

	# Filter by common markers, samples and remove duplicated and NA phenos
	common        = filterByCommonMarkersSamples (genotypeFile, phenotypeFile)
	genotypeFile  = common$genotypeFile
	phenotypeFile = common$phenotypeFile
	trait         = common$trait

	# Print info trait, N samples
	msgmsg("Evaluating following trait: ", trait) 
	nSamples = ncol (common$genotype)
	msgmsg ("N =",nSamples,"individuals with phenotypic and genotypic information \n")

	if (config$filtering == TRUE) {
		msgmsg ("Using filters")
		# Apply filters to genotype (markers and samples) by calling external program
		filtered      = filterByQCFilters (common$genotypeFile, common$phenotypeFile, config)
		genotypeFile  = filtered$genotypeFile
		phenotypeFile = filtered$phenotypeFile

		# Remove no polymorohic markers (MAF > Threshold) and Write chromosome info (map.tbl)
		maf           = filterByMAF (genotypeFile, config$ploidy, config$MAF)
		genotypeFile  = maf$genotypeFile
	}else {
		msgmsg ("Without filters")

		# Remove no polymorohic markers (MAF > 0.0) and Write chromosome info (map.tbl)
		maf           = filterByMAF (genotypeFile, config$ploidy, 0.0)
		genotypeFile  = maf$genotypeFile

		# Create GWASpoly files
		runCommand (sprintf ("ln -s %s %s", basename (genotypeFile), "out/filtered-gwasp4-genotype.tbl"))
		runCommand (sprintf ("ln -s %s %s", basename (phenotypeFile), "out/filtered-gwasp4-phenotype.tbl"))

		# Create plink files
		plinkFile  = paste0 ("out/", strsplit (basename(genotypeFile), split="[.]")[[1]][1], "-plink")
		gwaspToPlinkFormat (genotypeFile, plinkFile)

		# Make plink binary files from text file
		cmm = paste ("plink --file", plinkFile, "--allow-extra-chr --make-bed", "--out", plinkFile)
		runCommand (cmm, "log-filtering.log")

		# Make links of no-filtered to filtered files 
		runCommand (sprintf ("ln -s %s.ped out/filtered-plink-genotype.ped", basename (plinkFile)), "log-filtering.log")
		runCommand (sprintf ("ln -s %s.map out/filtered-plink-genotype.map", basename (plinkFile)), "log-filtering.log")

		runCommand (sprintf ("ln -s %s.bed out/filtered-plink-genotype.bed", basename (plinkFile)), "log-filtering.log")
		runCommand (sprintf ("ln -s %s.fam out/filtered-plink-genotype.fam", basename (plinkFile)), "log-filtering.log")
		runCommand (sprintf ("ln -s %s.bim out/filtered-plink-genotype.bim", basename (plinkFile)), "log-filtering.log")
	}

	msgmsg ("Converting phenotype for plink and TASSEL and create VCF genotype for TASSEL...")
	gwasp2plinkPhenotype  (phenotypeFile,"out/filtered-plink-phenotype.tbl") 
	gwaspToTasselPhenotype (phenotypeFile,"out/filtered-tassel-phenotype.tbl") 
	plinkToVCFFormat ("out/filtered-plink-genotype", "out/filtered-tassel-genotype")

	return (list (genotypeFile=genotypeFile, phenotypeFile=phenotypeFile, trait=common$trait))
}
#-------------------------------------------------------------
# Return the format type of genotype
# Checks if VCF, GWASpoly(k-matrix-chrom-pos), k-matrix, and fitPoly
#-------------------------------------------------------------
genotypeProcessing <- function (genotypeFile, type, mapFile, ploidy) {
	msg ("Processing genotype file format...")

	if (type=="gwaspoly"){
		msgmsg ("Using gwaspoly genotype format..")
		newGenotypeFile = genotypeFile

	}else if (type=="kmatrix") {# Only for tetraploids
		msgmsg ("Converting kmatrix gwaspoly genotype...")
		newGenotypeFile = createGwaspolyGenotype (genotypeFile, mapFile)

	}else if (type=="vcf") {
		msgmsg ("Converting VCF genotype ...")
		newGenotypeFile = convertVCFToACGTByNGSEP (genotypeFile) #output: filename.csv

	}else if (type=="fitpoly"){
		newGenotypeFile = convertFitpolyToGwaspolyGenotype (genotypeFile, mapFile) #output: filename.csv

	}else {
		msgmsg ("Error: Unknown genotype file format")
		return (NULL)
	}

	return (newGenotypeFile)
}


#-------------------------------------------------------------
# Filter by missing markers and samples, MAF, and HWE
# Apply filters to genotype (markers and samples) by calling external program
#-------------------------------------------------------------
filterByQCFilters <- function (genotypeFile, phenotypeFile, config) 
{
	# Format convertion from gwasp4 to plink2
	#msgmsg ("Converting gwaspoly to plink formats...")
	plinkFile  = paste0 ("out/", strsplit (basename(genotypeFile), split="[.]")[[1]][1], "-plink")
	gwaspToPlinkFormat (genotypeFile, plinkFile)

	#cmm = paste ("plink --file", plinkFile, "--make-bed", "--out", paste0(plinkFile,"-QC"))
	cmm = paste ("plink --file", plinkFile, "--make-bed --allow-extra-chr", "--out", paste0(plinkFile,"-QC"))

	# Create string for plink command with filters 
	msgmsg ("Filtering by missingness, MAF, and HWE")
	# Filter missingness per sample (MIND)"
	if (!is.null(config$MIND)) cmm=paste (cmm, paste ("--mind", config$MIND))
	# Filter missingness per SNP    (GENO)
	if (!is.null(config$GENO)) cmm=paste (cmm, paste ("--geno", config$GENO))

	# Obsolete, calculate directly 
	#### Filter SNPs with a low minor allele frequency (MAF)
	####if (!is.null(config$MAF)) cmm=paste (cmm, paste ("--maf", config$MAF))

	# Filter SNPs which are not in Hardy-Weinberg equilibrium (HWE).
	if (!is.null(config$HWE)) cmm=paste (cmm, paste ("--hwe", config$HWE))

	# Recode to plink format adjusted for tassel and plink
	runCommand (cmm, "log-filtering.log" )
	cmm = sprintf ("plink --bfile %s-QC --allow-extra-chr --recode tab --out %s-QC", plinkFile, plinkFile)
	runCommand (cmm, "log-filtering.log")

#	# Check for linkage disequilibrium
#	msgmsg ("    >>>> Checking form linkage disequilibrium...")
#	cmm = sprintf ("plink --bfile %s-QC --indep-pairwise 50 10 0.9 --out %s-LD", plinkFile, plinkFile)
#	runCommand (cmm, "log-filtering.log")
#	cmm = sprintf ("plink --bfile %s-QC --exclude %s-LD.prune.out --make-bed --out %s-LDB", plinkFile, plinkFile, plinkFile)
#	runCommand (cmm, "log-filtering.log")
#	cmm = sprintf ("plink --bfile %s-LDB --recode tab --out %s-LD", plinkFile, plinkFile)
#	runCommand (cmm, "log-filtering.log")

	# Copy links of filtered plink files to main dir
	runCommand (sprintf ("ln -s %s-QC.ped out/filtered-plink-genotype.ped", basename (plinkFile)), "log-filtering.log")
	runCommand (sprintf ("ln -s %s-QC.map out/filtered-plink-genotype.map", basename (plinkFile)), "log-filtering.log")
	# Same for binaries
	runCommand (sprintf ("ln -s %s-QC.bed out/filtered-plink-genotype.bed", basename (plinkFile)), "log-filtering.log")
	runCommand (sprintf ("ln -s %s-QC.fam out/filtered-plink-genotype.fam", basename (plinkFile)), "log-filtering.log")
	runCommand (sprintf ("ln -s %s-QC.bim out/filtered-plink-genotype.bim", basename (plinkFile)), "log-filtering.log")

	# Get final markers and individuals"
	msgmsg ("Writting filtered markers and individuals...")
	filteredSamples <- as.character (read.table ("out/filtered-plink-genotype.ped", check.names=F)[,2])
	filteredMarkers <- as.character (read.table ("out/filtered-plink-genotype.map", check.names=F)[,2])

	# Filter phenotype (.tbl) with pheno names from QC filters 
	phenoAll = read.csv (phenotypeFile, header=T, check.names=F)
	rownames (phenoAll) = phenoAll [,1]
	phenoFiltered  = phenoAll [filteredSamples,]
	trait  <- colnames (phenoAll)[2]

	# Filter genotype (.tbl) with geno names from QC filters 
	genoAll  = read.csv (genotypeFile, header=T, check.names=F)
	rownames (genoAll) = genoAll [,1]
	filteredColumns = c (colnames (genoAll)[1:3], filteredSamples)
	genoFiltered <- genoAll [filteredMarkers, filteredColumns]

	msgmsg("Writting geno/pheno filtered by MAF, Missing, HWE...")
	outGenoFile  <- "out/filtered-gwasp4-genotype.tbl"
	outPhenoFile <- "out/filtered-gwasp4-phenotype.tbl"
	write.table (file=outGenoFile, genoFiltered, row.names=F, quote=F, sep=",")
	write.table (file=outPhenoFile, phenoFiltered, row.names=F, quote=F, sep=",")

	return (list (genotypeFile=outGenoFile, phenotypeFile=outPhenoFile, trait=trait))
}

#-------------------------------------------------------------
# Get alternate allele from a row of alleles
#-------------------------------------------------------------
get.ref <- function(x) {
	y <- paste(na.omit(x),collapse="")
	ans <- apply(array(bases),1,function(z,y){length(grep(z,y,fixed=T))},y)
	if (sum(ans)>2) {stop("Error in genotype matrix: More than 2 alleles")}
	if (sum(ans)==2) {ref.alt <- bases[which(ans==1)]}
	if (sum(ans)==1) {ref.alt <- c(bases[which(ans==1)],NA)}
	
	return(ref.alt)
}
#-------------------------------------------------------------
# Impute NA alleles
#-------------------------------------------------------------
impute.mode <- function(x) {
	ix <- which(is.na(x))
	if (length(ix)>0) {
		x[ix] <- as.integer(names(which.max(table(x))))
	}
	return(x)
}
#-------------------------------------------------------------
# Calculate threshold to decide SNPs significance
#-------------------------------------------------------------
calculateThreshold <- function (level, scores, method="FDR") 
{
	scores <- as.vector(na.omit (scores))
	m <- length(scores)
	if (method=="Bonferroni") 
		threshold <- -log10(level/m)
	else if (method=="FDR") {
		tmp <- cbind(10^(-scores),.qvalue(10^(-scores)))
		tmp <- tmp[order(tmp[,2]),]
		if (tmp[1,2] > level) {
			threshold <- -log10(tmp[1,1])*1.2
		} else {
			k <- max(which(tmp[,2] < level))
			threshold <- -log10(mean(tmp[k:(k+1),1]))
		}
	}

	return (threshold)
}

.qvalue <- function(p) {
        smooth.df = 3
        if (min(p) < 0 || max(p) > 1) {
            print("ERROR: p-values not in valid range.")
            return(0)
        }
        lambda = seq(0, 0.9, 0.05)
        m <- length(p)
        pi0 <- rep(0, length(lambda))
        for (i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
        }

        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0 <- predict(spi0, x = max(lambda))$y
        pi0 <- min(pi0, 1)
        if (pi0 <= 0) {
            print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
            return(0)
        }
        u <- order(p)
        qvalue.rank <- function(x) {
            idx <- sort.list(x)
            fc <- factor(x)
            nl <- length(levels(fc))
            bin <- as.integer(fc)
            tbl <- tabulate(bin)
            cs <- cumsum(tbl)
            tbl <- rep(cs, tbl)
            tbl[idx] <- tbl
            return(tbl)
        }
        v <- qvalue.rank(p)
        qvalue <- pi0 * m * p/v
        qvalue[u[m]] <- min(qvalue[u[m]], 1)
        for (i in (m - 1):1) {
            qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 
                1)
        }
        return(qvalue)
    }
#-------------------------------------------------------------
# If genotype is VCF, convert to ACGT kmatrix (.csv)
#-------------------------------------------------------------
convertGenotypeVCFtoACGT <- function (genotypeFile) {
	msgmsg ("Checking genotype file format...")

	con = file(genotypeFile,"r")
	firstLine = readLines (con, n=1)
	close (con)

	if (grepl ("VCF", firstLine)) {
		msgmsg ("Converting VCF genotype to k-matrix genotype (.csv)")
		genotypeFile = convertVCFToACGTByNGSEP (genotypeFile) #output: filename.csv
	}
	return (genotypeFile)
}

#-------------------------------------------------------------
# Get common sample names
#-------------------------------------------------------------
filterByCommonMarkersSamples <- function (genotypeFile, phenotypeFile, mapFile=NULL) {
	msgmsg ("Filtering by common markers and samples...")
	geno  = read.csv (file=genotypeFile, check.names=F)
	pheno = read.csv (file=phenotypeFile, check.names=F)

	# From phenotype, remove duplicated and NA samples
	pheno  = pheno [!is.na(pheno[,2]),]
	pheno  = pheno [!duplicated (pheno[,1]),]

	# From genotype, remove duplicated markers
	geno  = geno [!duplicated (geno[,1]),]

	# Get common samples
	genoSamples   = colnames (geno[,-(1:3)])
	phenoSamples  = pheno [,1] 
	commonSamples = intersect (genoSamples, phenoSamples) 

	# Create new geno, pheno
	genoCommon  = geno  [,c (colnames(geno)[1:3], commonSamples)]
	phenoCommon = pheno [pheno[,1] %in% commonSamples,]
	trait       = colnames(phenoCommon)[2]

	genoCommonFile  = paste0 ("out/", addLabel (genotypeFile, "COMMON"))
	phenoCommonFile = paste0 ("out/", addLabel (phenotypeFile, "COMMON"))
	write.csv (file=genoCommonFile, genoCommon, quote=F, row.names=F)
	write.csv (file=phenoCommonFile, phenoCommon, quote=F, row.names=F)

	return (list (genotypeFile=genoCommonFile, phenotypeFile=phenoCommonFile, trait=trait))
}
#-------------------------------------------------------------
# Impute, filter by MAF, unify geno and pheno names
# Only for "ACGT" format (For other formats see GWASpoly sources)
#-------------------------------------------------------------
filterByMAF <- function(genotypeFile, ploidy, thresholdMAF){
	msgmsg ("Reading phenotype for MAF processing....")
	if (is.null (thresholdMAF)) thresholdMAF = 0.0

	geno              <- read.csv(genotypeFile,check.names=F, as.is=T)
	gid.geno          <- colnames(geno)[-(1:3)]
	markers           <- as.matrix(geno[,-(1:3)])
	rownames(markers) <- geno[,1]
	
	tmp     <- apply(markers,1,getReferenceAllele)
	map     <- data.frame(Marker=geno[,1],Chrom=factor(geno[,2],ordered=T),Position=geno[,3],stringsAsFactors=F)
	map$Ref <- tmp[1,]
	map$Alt <- tmp[2,]

	msgmsg ("Calculating numeric genotype matrix...")
	matRefMarkers = cbind (map$Ref, markers)

	#-------------------------------
	# Convert All ACGT matrix to Num
	acgtToNum <- function(x){
		y <- gregexpr(pattern=x[1],text=x[-1],fixed=T)  
		ans <- as.integer(lapply(y,function(z){ifelse(z[1]<0,ploidy,ploidy-length(z))}))	
		return(ans)
	}

	matTransposed   = t(matRefMarkers)
	ACGTList        = mclapply(seq_len(ncol(matTransposed)), function(i) matTransposed[,i],mc.cores=NCORES)
	numList         = mclapply(ACGTList, acgtToNum, mc.cores=NCORES)
	M               = as.data.frame (numList,col.names=rownames (matRefMarkers))
	rownames(M)     = gid.geno
	#-------------------------------

	# Check LG Global MAF (AF?)
	msgmsg ("Checking minor allele frecuency, MAF=", thresholdMAF)
	MAF         <- apply(M,2,function(x) {AF <- mean(x,na.rm=T)/ploidy;MAF <- ifelse(AF > 0.5,1-AF,AF)})
	polymorphic <- which(MAF>thresholdMAF)

	M   <- M[,polymorphic]
	map <- map[polymorphic,]
	map <- map[order(map$Chrom,map$Position),]
	M   <- M[,map$Marker]
	m   <- nrow(map)
	msgmsg ("Number of polymorphic markers:",m,"\n")
	
	missing <- which(is.na(M))
	if (length(missing)>0) {
		msgmsg("Missing marker data imputed with population mode...")
		M <- apply(M,2,impute.mode)
	}

	# Write geno MAF
	rownames (geno) = geno [,1]
	genoMAF         = geno [colnames(M),]
	genoMAFFile    = paste0 ("out/", addLabel (basename(genotypeFile), "MAF"))
	write.csv (genoMAF, genoMAFFile, quote=F, row.names=F)

	# Write chromosome info 
	write.table (file="out/map.tbl", map, quote=F, row.names=F, sep="\t")

	return (list (genotypeFile=genoMAFFile, geno=genoMAF))
}

#-------------------------------------------------------------
# Impute, filter by MAF, unify geno and pheno names
# Only for "ACGT" format (For other formats see GWASpoly sources)
#-------------------------------------------------------------
old_filterByMAFCommonNames <- function(genotypeFile, phenotypeFile, ploidy){
	NCORES    = detectCores()
	format    = "ACGT"	
	n.traits  = 1
	bases     = c("A","C","G","T")

	thresholdMAF=0.0
	
	msgmsg ("Reading genotype and phenotype....")
	geno      <- read.csv(file=genotypeFile,as.is=T,check.names=F)
	gid.geno  <- colnames(geno)[-(1:3)]
	pheno     <- na.omit (read.csv(file=phenotypeFile,as.is=T,check.names=F))
	gid.pheno <- unique(pheno[,1])

	msgmsg ("Removing duplicated markers...")
	geno <- geno [!duplicated (geno[,1]),]  
	msgmsg ("Removing duplicated samples...")
	pheno <- pheno [!duplicated (pheno[,1]),]  

	markers <- as.matrix(geno[,-(1:3)])
	rownames(markers) <- geno[,1]
	
	tmp     <- apply(markers,1,getReferenceAllele)
	map     <- data.frame(Marker=geno[,1],Chrom=factor(geno[,2],ordered=T),Position=geno[,3],stringsAsFactors=F)
	map$Ref <- tmp[1,]
	map$Alt <- tmp[2,]

	msgmsg ("Calculating numeric genotype matrix...")
	matRef = cbind (map$Ref, markers)

	#-------------------------------
	# Convert All ACGT matrix to Num
	acgtToNum <- function(x){
		y <- gregexpr(pattern=x[1],text=x[-1],fixed=T)  
		ans <- as.integer(lapply(y,function(z){ifelse(z[1]<0,ploidy,ploidy-length(z))}))	
		return(ans)
	}

	matRefMarkers   = matRef
	matTransposed   = t(matRefMarkers)
	ACGTList        = mclapply(seq_len(ncol(matTransposed)), function(i) matTransposed[,i],mc.cores=NCORES)
	numList         = mclapply(ACGTList, acgtToNum, mc.cores=NCORES)
	M               = as.data.frame (numList,col.names=rownames (matRefMarkers))
	rownames(M)     = gid.geno
	#-------------------------------

	# Check LG Global MAF (AF?)
	msgmsg ("Checking minor allele frecuency, MAF=", thresholdMAF)
	MAF <- apply(M,2,function(x){AF <- mean(x,na.rm=T)/ploidy;MAF <- ifelse(AF > 0.5,1-AF,AF)})
	polymorphic <- which(MAF>thresholdMAF)

	M <- M[,polymorphic]
	map <- map[polymorphic,]
	map <- map[order(map$Chrom,map$Position),]
	M <- M[,map$Marker]
	m <- nrow(map)
	msgmsg("Number of polymorphic markers:",m,"\n")
	
	missing <- which(is.na(M))
	if (length(missing)>0) {
		msgmsg("Missing marker data imputed with population mode...")
		M <- apply(M,2,impute.mode)
	}
	
	msgmsg("Matching genotypic and phenotypic data...")
	gid <- intersect(gid.pheno, gid.geno)
	pheno <- pheno[is.element(pheno[,1],gid),]
	rownames (pheno) = pheno [,1]
	M <- M[gid,]
	N <- length(gid)
	msgmsg ("N =",N,"individuals with phenotypic and genotypic information \n")
	
	n.fixed <- ncol(pheno) - n.traits - 1
	if (n.fixed > 0) {
		fixed <- data.frame(pheno[,(n.traits+2):ncol(pheno)],stringsAsFactors=F)
		fixed.names <- colnames(pheno)[(n.traits+2):ncol(pheno)]
		colnames(fixed) <- fixed.names
		pheno <- data.frame(pheno[,1:(1+n.traits)],stringsAsFactors=F)
		cat(paste("Detected following fixed effects:\n",paste(fixed.names,collapse="\n"),"\n",sep=""))
	} else {
		fixed <- data.frame(NULL)
	}
	trait <- colnames(pheno)[-1]
	msgmsg("Evaluating following trait: ", trait) 
	

	#msgmsg ("Writing geno/pheno filtered by MAF, duplicated, common names")
	rownames (geno) = geno [,1]
	genoCommon      = geno [colnames(M),c(colnames(geno)[1:3], pheno[,1])]
	phenoCommon     = pheno 
	genoCommonFile  = paste0 ("out/", addLabel (genotypeFile, "COMMON"))
	phenoCommonFile = paste0 ("out/", addLabel (phenotypeFile, "COMMON"))
	write.csv (file=genoCommonFile, genoCommon, quote=F, row.names=F)
	write.csv (file=phenoCommonFile, phenoCommon, quote=F, row.names=F)
	# Write chromosome info 
	map = genoCommon [, (1:3)]
	write.table (file="out/map.tbl", map, quote=F, row.names=F, sep="\t")
	# construct GWASpoly data structure
	gwaspolyData = new("GWASpoly",map=map,pheno=pheno,fixed=fixed,geno=M,ploidy=ploidy)

	return (list (genotypeFile=genoCommonFile, phenotypeFile=phenoCommonFile, data=gwaspolyData, trait=trait))
}

#-------------------------------------------------------------
# Add label to filename and new extension (optional)
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
#-------------------------------------------------------------
# Add label to filename
#-------------------------------------------------------------
old_addLabel <- function (filename, label)  {
	nameext = strsplit (filename, split="[.]")
	newName = paste0 (nameext [[1]][1], "-", label, ".", nameext [[1]][2])
	return (newName)
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
  messages = unlist (list (...))
  cat ("\n>>>>", messages, "\n")
}

msgmsg <- function (...) 
{
  messages = unlist (list (...))
  cat ("\t>>>>", messages, "\n")
}

msgError <- function (...) {
		messages = unlist (list (...))
		cat (strrep("-", sum(sapply(messages, nchar))),"\n")
		cat (messages, "\n")
		cat (strrep("-", sum(sapply(messages, nchar))),"\n")
}

#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
view <- function (data, n=5,m=6) {
	name = paste (deparse (substitute (data)),":  ")
	if (is.null (dim (data))) {
		dimensions = paste (length (data))
		message (name, "(", paste0 (dimensions),")")
		if (length (data) < 6) n = length(data)
		print (data[1:n])
	}else {
		dimensions = paste0 (unlist (dim (data)),sep=c(" x ",""))
		message (name, "(", paste0 (dimensions),")")
		if (nrow (data) < 5) n = nrow(data)
		if (ncol (data) < 6) m = ncol(data)
		print (data[1:n,1:m])
	}
	#write.csv (data, paste0("x-", filename, ".csv"), quote=F, row.names=F)
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
# Call main 
#-------------------------------------------------------------
withCallingHandlers (
	main (), 
	error = function (e) { 
		if (DEBUG){
			print (sys.calls()[-1])
			quit ()
		}
	},
	warning = function (w) { 
		message (geterrmessage ())
	}
)
