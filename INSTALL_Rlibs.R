#
# Create multiGWAS profile according to current directory
# 

# Create multiGWAS profile
message ("Creating multiGWAS profile...")

sink ("multiGWAS_profile.sh", append=F)
writeLines ("")
writeLines ("#------------------- multiGWAS.R profile ---------------------")
MULTIGWAS_HOME = getwd ()
#writeLines ('echo "...multiGWAS profile"')
writeLines (paste0 ("export MULTIGWAS_HOME=", getwd()))
writeLines ("MULTIGWAS_TOOLS=$MULTIGWAS_HOME/tools")
writeLines ("MULTIGWAS_SOURCES=$MULTIGWAS_HOME/sources")
writeLines ("export PATH=$PATH:$MULTIGWAS_TOOLS:$MULTIGWAS_SOURCES")
sink()


# Write into .bashrc
profileFile = paste0 (path.expand ("~"), "/.bashrc")
sink (profileFile, append=T)
writeLines ("")
writeLines ("#------------------- multiGWAS.R tool profile ---------------------")
writeLines (paste0 (". ", getwd(), "/multiGWAS_profile.sh"))
sink ()

message ("")
message ("MultiGWAS is ready to use, right after installed!")
message ("")

# Install R libraries
libpath = paste0 (MULTIGWAS_HOME, "/opt/Rlibs")
message (libpath)
message ("\n\nInstalling R libraries...\n\n")

.libPaths (libpath)

if (!require("pacman")) install.packages('pacman', lib=libpath, repos='http://cran.us.r-project.org')

pacman::p_load("rrBLUP", "parallel","config","dplyr","stringi","qqman","VennDiagram","RColorBrewer","circlize","gplots", "rmarkdown", "kableExtra") 

if (!require("GWASpoly")) install.packages('opt/GWASpoly_1.3.tar.gz', lib=libpath, repos=NULL, type="source")

