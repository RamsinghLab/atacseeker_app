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
install.packages(c("devtools", "knitr", "locfit", "preseqR", 
    "R2HTML", "rjson", "rmarkdown", "statmod"))

## Bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite(pkgs=c("BSgenome.Hsapiens.UCSC.hg19", "csaw", "edgeR", "Gviz", "Homo.sapiens", 
    "LOLA", "PWMEnrich", "PWMEnrich.Hsapiens.background", "Rsamtools", "rtracklayer"))

## Github package
library(devtools)
install_github("RamsinghLab/ATACseeker")
