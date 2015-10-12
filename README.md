# ATACseeker

Repo for BaseSpace ATACseq pipeline. 

*requires*: docker 
*requires*: wget, git, pandoc, [samblaster](https://github.com/GregoryFaust/samblaster), samtools, [bwa](https://github.com/lh3/bwa)
*requires*: R, knitr, [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html), statmod, locfit, edgeR, [csaw](http://bioconductor.org/packages/release/bioc/html/csaw.html)

This repository contains a first pass attempt for an application for BaseSpace.

The elements of the approach are as follows: 

* Docker image with all packages (*bwa, samtools, samblaster etc.*) installed. 
* R wrapper for the alignment pipeline. 
* Quality checks on aligned reads.
* Pipe reads through `csaw` for differential accessibility testing.
* Deploy as Native App on BaseSpace. 

## Directories

Main development directory for BaseSpace App is `nativeApp`.

`csaw`: Testing code for csaw pipeline.
`docker`: Testing code for Docker Image. This docker image only does the alignment of raw reads.
`scripts`: random scripts live here.
