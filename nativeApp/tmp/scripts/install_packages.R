##
# Set CRAN mirror and install R packages
##

## Hard code the US repo for CRAN
r <- getOption("repos")             
r["CRAN"] <- "http://cran.us.r-project.org"
options(repos = r)
rm(r)

##
# Install required packages
##

## CRAN packages
## install.packages('knitr')
## install.packages('rmarkdown')
## install.packages('rjson')
## install.packages('statmod')
## install.packages('locfit')

install.packages(c("knitr", "rmarkdown", "rjson", "statmod", "locfit"))

## Bioconductor packages
source("https://bioconductor.org/biocLite.R")
## biocLite('edgeR')
## biocLite('csaw')
## biocLite('Homo.sapiens')
## biocLite('Gviz')
biocLite(pkgs=c("Rsamtools", "rtracklayer", "edgeR", "csaw", "Homo.sapiens", "LOLA"))
