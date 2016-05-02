#/bin/bash

BASE=".."

if [ $# -eq 0 ]; then
    echo "Usage: bash $0 <in.bam>"
    exit 0
fi

chrM_bam=$1
in_base=`basename ${chrM_bam} .bam`
dir=`dirname ${chrM_bam}`

samtools index ${chrM_bam}
samtools sort ${chrM_bam} ${dir}/${in_base}.sorted
samtools rmdup ${dir}/${in_base}.sorted.bam ${dir}/${in_base}.sorted.nodup.bam
samtools index ${dir}/${in_base}.sorted.nodup.bam

python assembleMTgenome.py \
-f ${BASE}/reference/chrM.fa \
-i ${dir}/${in_base}.sorted.nodup.bam \
-a ${BASE}/reference/hg19.fa \
-q 25 -c 0.75 -d 5 -s `which samtools` -FCUP \
-o ${dir}/${in_base}
