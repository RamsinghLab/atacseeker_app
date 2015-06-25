## check quality of reads in fastq files

fastqc normal1_1.1mReads.fastq
fastqc normal1_2.1mReads.fastq

## index the reference
bwa index hg19.fa

## align paired end reads by using bwa mem and pipe it to samblaster to remove duplicates. Use samtools to make BAM file 
bwa mem -M hg19.fa normal1_2.1mReads.fastq normal1_2.1mReads.fastq | samblaster | samtools view -Sb > normal1.bam

## sort and index BAM for downstream processing
samtools sort normal1.bam normal1_sorted
samtools index normal1_sorted

## filter out MT DNA
samtools idxstats normal1_sorted.bam | cut -f 1 | grep -v chrM | xargs samtools view -b normal1_sorted.bam > normal1_sorted_filtered.bam