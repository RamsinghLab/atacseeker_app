# ATACseq_counts
Repo for BaseSpace ATACseq pipeline. 

This repository contains a first pass attempt for an application for BaseSpace.

The elements of the approach are as follows: 

1. Write an R wrapper for the alignment pipeline. 
2. Create a docker image with all packages (*bwa, samtools, samblaster etc.*) installed. 
3. Deploy this as a Native App on BaseSpace. 

The long term objective is to integrate this with the **ATACseeker** package and **csaw**.

*This README will be edited over the course of development.*
