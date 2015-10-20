##
# Set CRAN mirror
##

r <- getOption("repos")             # hard code the US repo for CRAN
r["CRAN"] <- "http://cran.us.r-project.org"
options(repos = r)
rm(r)

##
# Install required packages
##

install.packages('knitr')
install.packages('rmarkdown')
install.packages('rjson')
install.packages('statmod')
install.packages('locfit')
install.packages('edgeR')
install.packages('csaw')


