# ATACseeker #

Illumina BaseSpace ATACseq pipeline. 

This repository contains code for a Illumina BaseSpace Native App for doing ATACseq analysis.

The elements of the approach are as follows: 

* Docker image with all packages (*bwa, samtools, samblaster etc.*) installed. 
* R wrapper for the alignment pipeline. 
* Quality checks on aligned reads.
* Pipe reads through `csaw` for differential accessibility testing.
* Deploy as Native App on BaseSpace. 

## Directories

Main development directory for BaseSpace App is `nativeApp`.

* `csaw`: Testing code for csaw pipeline.
* `docker`: Testing code for docker image. This image only implements the alignment pipeline.
* `scripts`: random scripts live here.

*requires*: 
* docker, wget, pandoc
* R
  * CRAN: knitr, [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html), rjson, statmod, locfit
  * Bioconductor: edgeR, [csaw](http://bioconductor.org/packages/release/bioc/html/csaw.html)
