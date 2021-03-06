#!/bin/bash -x

if [ $# -eq 0 ]; then
	echo "Usage: bash $0 reference fastq_prefix [genome]"
	exit 0
fi

GENOME=$1
SAMPLE=$2
SUFFIX="fq"
FASTQS=$SAMPLE".R1."$SUFFIX" "$SAMPLE".R2."$SUFFIX 
GENOME_TYPE=${3:-"chrM"} ## let's hope we know this in advance!!
# CHRSIZES=$GENOME".chrom.sizes" ## assuming we do

##
## align paired end reads by using bwa mem,
## pipe it to samblaster to remove duplicates, 
## and use samtools to make BAM file
## 
## flags: no secondary, no QC failures, minimum quality 10, 
##        SAM format in, BAM format out
## 

FLAGS="-F 0x100i -F 0x200i -F 0x400i -F 0x800i -q 10 -Sb"

##
## spell out these conditions in the resulting BAM file name 
## 

OUT_BAM=$SAMPLE".unique_q10minimum."$GENOME_TYPE".bam"

##
## do the alignment
##

echo "bwa mem  $GENOME $FASTQS | \
  samblaster  | \
  samtools view $FLAGS - > $OUT_BAM"

bwa mem $GENOME $FASTQS | \
  samblaster | \
  samtools view $FLAGS - > $OUT_BAM

##
## sort and index the BAM for downstream processing
##

SORTED=`basename $OUT_BAM .bam`.sorted.bam
DIR=`dirname $SAMPLE`
samtools sort $OUT_BAM -o $DIR"/"$SORTED
samtools index $DIR"/"${SORTED}

echo "alignment file created ..." ${SORTED}
