# ATACseeker

Repo for BaseSpace ATACseq pipeline. 

*requires*: R, bwa, samtools, samblaster, genomeCoverageBed, wget, mysql

This repository contains a first pass attempt for an application for BaseSpace.

The elements of the approach are as follows: 

* Docker image with all packages (*bwa, samtools, samblaster etc.*) installed. 
* R wrapper for the alignment pipeline. 
* Quality checks on aligned reads.
* Pipe reads through *csaw* for differential accessibility testing.
* Deploy as Native App on BaseSpace. 
