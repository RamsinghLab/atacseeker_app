#!/bin/bash

if [ $# -eq 0 ]; then
	echo "Usage: bash $0 <in.bam> dedup?"
	exit 0
fi

in_bam=$1
in_base=`basename ${in_bam} .bam`
dir=`dirname ${in_bam}`
OUT_DIR=${3:-"/data/scratch"}

if [[ -n $2 ]]; then
	samtools rmdup ${dir}/${in_base}.bam ${OUT_DIR}/${in_base}.nodup.bam
	samtools index ${OUT_DIR}/${in_base}.nodup.bam
	in_bam="${OUT_DIR}/${in_base}.nodup.bam"
	in_base=`basename ${in_bam} .bam`
	dir="${OUT_DIR}"
fi


GENOME=${4:-"hg19"}
CHROM_SIZES="chrom.sizes/"$GENOME".chrom.sizes"

## if we don't already have chrom.sizes for the genome in question, fetch them
## test -f $CHROM_SIZES || ./fetchChromSizes $GENOME > $CHROM_SIZES
echo -n "Creating coverage bedgraph for ${in_base} ... "
genomeCoverageBed -bg -split -ibam ${dir}/${in_base}.bam -g $CHROMSIZES > ${OUT_DIR}/${in_base}.cov.bg && \
	sort -k1,1 -k2,2n ${OUT_DIR}/${in_base}.cov.bg > ${OUT_DIR}/${in_base}.cov.sorted.bg && \
	bedGraphToBigWig ${OUT_DIR}/${in_base}.cov.sorted.bg $CHROM_SIZES ${OUT_DIR}/${in_base}.bw && \
	rm ${OUT_DIR}/${in_base}.cov.bg  ${OUT_DIR}/${in_base}.cov.sorted.bg
echo "...done."
