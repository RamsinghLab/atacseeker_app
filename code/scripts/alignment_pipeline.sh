## > For expository purposes, let's generalize a tiny bit:
##
SUBJECT="normal1"
SUFFIX="1mReads.fastq"
FASTQS=$SUBJECT"_1."$SUFFIX" "$SUBJECT"_2."$SUFFIX 
## echo $FASTQS
## > normal1_1.1mReads.fastq normal1_2.1mReads.fastq
 
GENOME="hg19" ## let's hope we know this in advance!!
CHRSIZES=$GENOME".chrom.sizes" ## assuming we do

## check quality of reads in fastq files
## 
## > good idea! but you can do this separately in BaseSpace with its own app --t
## > if we were going to do it locally, a loop would make life a little easier:
## 
## for FASTQ in $FASTQS; do
##   fastqc $FASTQ
## done

## index the reference
##
## > also a good idea, but should be done beforehand & uploaded with the app --t
##
## bwa index $GENOME.fa

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
bwa mem -M $GENOME.fa $FASTQS | samblaster | grep -v "chrM" | \
  samtools view $FLAGS - >

## sort and index the BAM for downstream processing
##
SORTED=`basename $OUTBAM .bam`.sorted
samtools sort $OUTBAM $SORTED 
samtools index $SORTED".bam"

## if you wanted to, you could bamToBigWig it at this point:
##
scripts/bamToBigWig $SORTED $GENOME
