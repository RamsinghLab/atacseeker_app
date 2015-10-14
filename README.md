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
* `docker`: Testing code for docker image. This docker image only implements the alignment pipeline.
* `scripts`: random scripts live here.

*requires*: 
* docker 
* wget, git, pandoc, [samblaster](https://github.com/GregoryFaust/samblaster), samtools, [bwa](https://github.com/lh3/bwa)
* R, knitr, [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html), statmod, locfit, edgeR, [csaw](http://bioconductor.org/packages/release/bioc/html/csaw.html)
