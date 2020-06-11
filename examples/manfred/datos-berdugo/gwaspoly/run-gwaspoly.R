#!/usr/bin/Rscript
library (GWASpoly)

# Run GWASpoly (Naive and Full) using as input
# a genotype (k-matrix) and phenotype

# Paramters 
TRAIT     = "train_CPPT"
N_TRAITS  = 1
FORMAT    = "ACGT"   # "AB" "ACGT" "numeric"
STRUCTURE = T
NCORES    = 8

#-------------------------------------------------------------
# Add label to filename
#-------------------------------------------------------------
addLabel <- function (filename, label)  {
	nameext = strsplit (filename, split="[.]")
	newName = paste0 (nameext [[1]][1], "-", label, ".", nameext [[1]][2])
	return (newName)
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) {
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}
#-------------------------------------------------------------
  
args = commandArgs(trailingOnly = TRUE)
print (args)

genotypeFile  = args [1]
phenotypeFile = args [2]

#genotypeFile  = "agrosavia-genotype-ACGT-CLEANED-MAP-COMMON.csv"
#phenotypeFile = "agrosavia-phenotype-GOTA-CLEANED-COMMON.csv"


message (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
message ("Genotype file  : ", genotypeFile)
message ("Phenotype file : ", phenotypeFile)
message ("Trait          : ", TRAIT) 
message ("Format         : ", FORMAT) 
message ("Structure      : ", STRUCTURE) 
message (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

ph = read.table (phenotypeFile, header=T, sep=",")
gn = read.table (genotypeFile, header=T, sep=",")


# Read input genotype and genotype (format: "numeric" or "ACGT")
message (">>> Reading data...")
data = read.GWASpoly (ploidy = 4, delim=",", format = FORMAT, n.traits = N_TRAITS, 
                      pheno.file = phenotypeFile, geno.file = genotypeFile)

# Populations structure by kinship
# Used to include population structure covariates
if (STRUCTURE==T) {
	message (">>> Full: With structure")
	PARAMS <- set.params(n.PC=10)
	#PARAMS <- set.params(fixed=c("K1","K2","K3","K4","K5"), fixed.type=rep("numeric",5))
	data2 <- set.K(data)
} else {
	message (">>> Naive: No structure")
	PARAMS = NULL
	data2 <- set.K(data, K=NULL)
}

#PARAMS <- set.params(fixed=c("K1","K2","K3","K4"),
#                     fixed.type=rep("numeric",4))
# GWAS execution
message (">>> Running GWASpoly...")
data3 = GWASpoly(data2, models=c("general","additive","1-dom", "2-dom"),
				 traits=c(TRAIT),n.core=NCORES, params=PARAMS)

# QTL Detection
data4 = set.threshold (data3, method="Bonferroni",level=0.05, n.core=NCORES)
#data4 = set.threshold (data3, method="FDR",level=0.05, n.core=4)
significativeQTLs = get.QTL (data4)
outFile = paste0 ("out-", TRAIT, "-QTL.csv")
write.csv (significativeQTLs, outFile)

# Write all results
outFile = paste0 ("out-", TRAIT, "-results.csv")
write.GWASpoly(data4, TRAIT, outFile, what = "scores", delim = ",")

# Plots
models <- c("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
message (">>> Ploting results...")
outFile = paste0 ("out-manhattan-QQ-", TRAIT,".pdf")

# Manhattan plot Output
pdf (file=outFile, width=7, height=14)
  par(mfrow=c(4,3)) #specifies a 1 x 3 panel
  for (i in 1:6) {
	MODEL = models [i]
	message (">>> Manhattan Plot for model: ", MODEL)
    manhattan.plot (data4, trait=TRAIT, model=MODEL)
  }  
#dev.off ()

# QQ-plot Output
outFile = paste0 ("out.", TRAIT,"-QQ.pdf")
#pdf (file=outFile)
  #par(mfrow=c(2,3)) #specifies a 2 x 3 panel
  for (i in 1:6) {
	  MODEL = models [i]
	  message (">>> QQ Plot for model: ", MODEL)
    qq.plot(data3,trait=TRAIT,model=models[i])
  }
dev.off()



