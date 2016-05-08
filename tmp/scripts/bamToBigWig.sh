#!/bin/bash

if [ $# -eq 0 ]; then
	echo "Usage: bash $0 <in.bam>"
	exit 0
fi

in_bam=$1

in_base=`basename ${in_bam} .bam`
dir=`dirname ${in_bam}`
GENOME=${2:-"hg19"}
CHROM_SIZES="chrom.sizes/"$GENOME".chrom.sizes"

## if we don't already have chrom.sizes for the genome in question, fetch them
## test -f $CHROM_SIZES || ./fetchChromSizes $GENOME > $CHROM_SIZES

echo -n "Creating coverage bedgraph for ${in_base} ... "
genomeCoverageBed -bg -split -ibam ${dir}/${in_base}.bam -g $CHROMSIZES > ${dir}/${in_base}.cov.bg && \
	sort -k1,1 -k2,2n ${dir}/${in_base}.cov.bg > ${dir}/${in_base}.cov.sorted.bg && \
	./bedGraphToBigWig ${dir}/${in_base}.cov.sorted.bg $CHROM_SIZES ${dir}/${in_base}.bw && \
	rm ${dir}/${in_base}.cov.bg  ${dir}/${in_base}.cov.sorted.bg
echo "...done."
