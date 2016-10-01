#/bin/bash -xv

if [ $# -eq 0 ]; then
    echo "Usage: bash $0 mt_dir reference"
    exit 0
fi 

mt_dir=$1
ref=$2

bamaddrg -c $(for b in ${mt_dir}/*/*.bam; do echo -b "$b" -s $(basename "$b" .unique_q10minimum.chrM.sorted.bam) -r $(basename "$b" .unique_q10minimum.chrM.sorted.bam); done) > ${mt_dir}/all_bams.bam
samtools index ${mt_dir}/all_bams.bam
freebayes -C 5 -p 1 -f ${ref} ${mt_dir}/all_bams.bam | vcffilter -f "QUAL >20" > ${mt_dir}/var.vcf
cat ${mt_dir}/var.vcf | vcf-to-tab > ${mt_dir}/var.tsv
rm -rf ${mt_dir}/*.bam ${mt_dir}/*.bai
