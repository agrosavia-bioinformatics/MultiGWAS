#
# Create multiGWAS profile according to current directory
# 

MULTIGWAS_HOME = strsplit (getwd (), "install")[[1]][1]
message ("MULTIGWAS_HOME=",MULTIGWAS_HOME)

message ("Creating multiGWAS profile...")
sink ("multiGWAS_profile.sh", append=F)
writeLines ("\n#------------------- multiGWAS.R profile ---------------------")
writeLines (paste0 ("export MULTIGWAS_HOME=", MULTIGWAS_HOME))
writeLines ("MULTIGWAS_TOOLS=$MULTIGWAS_HOME/tools")
writeLines ("MULTIGWAS_SOURCES=$MULTIGWAS_HOME/sources")
writeLines ("export PATH=$PATH:$MULTIGWAS_TOOLS:$MULTIGWAS_SOURCES")
sink()


# Write into .bashrc
profileFile = paste0 (path.expand ("~"), "/.bashrc")
sink (profileFile, append=T)
writeLines ("")
writeLines ("#------------------- multiGWAS.R tool profile ---------------------")
writeLines (paste0 (". ", MULTIGWAS_HOME, "/multiGWAS_profile.sh"))
sink ()

message ("\nMultiGWAS is ready to use, right after installed!\n")

# Install R libraries
dir.create (paste0 (MULTIGWAS_HOME, "/opt/Rlibs"))
libpath = paste0 (MULTIGWAS_HOME, "/opt/Rlibs")
message (libpath)
message ("\n\nInstalling R libraries...\n\n")

.libPaths (libpath)

if (!require("pacman")) install.packages('pacman', lib=libpath, repos='http://cran.us.r-project.org')

pacman::p_load("rrBLUP", "parallel","config","dplyr","stringi","qqman","VennDiagram","RColorBrewer","circlize","gplots", "rmarkdown", "kableExtra" ,"doParallel", "ldsep") 

install.packages('opt/GWASpoly_1.3.tar.gz', lib=libpath, repos=NULL, type="source")

