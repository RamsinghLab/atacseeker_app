#!/bin/bash -xv

Rscript -e 'rmarkdown::render("/atacseeker/scripts/atacseeker.Rmd")'
cp /atacseeker/scripts/atacseeker.html '/data/output/appresults/27393372/ATACseeker 02-06-2016 11:17:27'
