#!/bin/bash -xv

cd /atacseeker/scripts
Rscript -e 'rmarkdown::render("atacseeker.Rmd")'
