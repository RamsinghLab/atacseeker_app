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
install.packages(c("knitr", "rmarkdown", "rjson", "statmod", "locfit"))

## Bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite(pkgs=c("Rsamtools", "rtracklayer", "edgeR", "csaw", "Homo.sapiens", "LOLA", "PWMEnrich", "PWMEnrich.Hsapiens.background", "BSgenome.Hsapiens.UCSC.hg19"))
