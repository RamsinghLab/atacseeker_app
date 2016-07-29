#!/bin/bash

if [ $# -eq 0 ]; then
	echo "Usage: bash $0 <in.bam> dedup?"
	exit 0
fi

in_bam=$1
in_base=`basename ${in_bam} .bam`
dir=`dirname ${in_bam}`

if [[ -z $2 ]]; then
	samtools rmdup ${dir}/${in_base}.bam ${dir}/${in_base}.nodup.bam
	samtools index ${dir}/${in_base}.nodup.bam
	in_bam="${dir}/${in_base}.nodup.bam"
	in_base=`basename ${in_bam} .bam`
fi

GENOME=${2:-"hg19"}
CHROM_SIZES="/atacseeker/scripts/chrom.sizes/"$GENOME".chrom.sizes"

## if we don't already have chrom.sizes for the genome in question, fetch them
## test -f $CHROM_SIZES || ./fetchChromSizes $GENOME > $CHROM_SIZES
echo -n "Creating coverage bedgraph for ${in_base} ... "
genomeCoverageBed -bg -split -ibam ${dir}/${in_base}.bam -g $CHROMSIZES > ${dir}/${in_base}.cov.bg && \
	sort -k1,1 -k2,2n ${dir}/${in_base}.cov.bg > ${dir}/${in_base}.cov.sorted.bg && \
	bedGraphToBigWig ${dir}/${in_base}.cov.sorted.bg $CHROM_SIZES ${dir}/${in_base}.bw && \
	rm ${dir}/${in_base}.cov.bg  ${dir}/${in_base}.cov.sorted.bg
echo "...done."
