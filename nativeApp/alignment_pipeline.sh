#!/bin/bash -x

DATA=$1

SUBJECT="normal1"
SUFFIX="1mReads.fastq"
FASTQS=$SUBJECT"_1."$SUFFIX" "$SUBJECT"_2."$SUFFIX 


REF_PATH="/atacseeker/reference/" 
GENOME="hg19" ## let's hope we know this in advance!!
CHRSIZES=$GENOME".chrom.sizes" ## assuming we do

##
## align paired end reads by using bwa mem,
## pipe it to samblaster to remove duplicates, 
## and use samtools to make BAM file
## 
## flags: no secondary, no QC failures, no PCR dupes, minimum quality 10, 
##        SAM format in, BAM format out
## 

FLAGS="-F 0x100i -F 0x200i -F 0x400i -F 0x800i -q 10 -Sb"

##
## spell out these conditions in the resulting BAM file name 
## 

OUTBAM=$SUBJECT".unique_q10minimum_noChrM."$GENOME".bam"

##
## a truly filthy trick: directly grep out the reads on chrM 
##

bwa mem -M $REF_PATH$GENOME.fa $FASTQS | \
  samblaster | \
  grep -v "chrM" | \
  samtools view $FLAGS - > $OUTBAM

##
## sort and index the BAM for downstream processing
##

SORTED=`basename $OUTBAM .bam`.sorted
samtools sort $OUTBAM $SORTED 
samtools index $SORTED".bam"

