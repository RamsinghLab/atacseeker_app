#!/bin/bash 

if [[ $# -eq 0 ]]; then
    echo "Usage: bash $0 <in.bam> output_dir reference"
    exit 0
fi

bam_file=$1
output_dir=$2
reference=$3

dir_name=`dirname ${bam_file}`
bam_name=`basename ${bam_file} .bam`
work_dir=${output_dir}/${bam_name}
mkdir -p ${work_dir}

externaltoolsfolder="/home/cmb-07/sn1/asifzuba/mtHeteroplasmy/ext_tools"

echo "**************************"
echo "Working with sample ..." ${bam_name}
echo "**************************"


echo "**************************"
echo "Extracting MT reads ..."
echo "**************************"

samtools view -F 1804 -q 10 -b ${bam_file} chrMT chrM > ${work_dir}/${bam_name}.MT.bam

echo "Done."
echo


echo "**************************"
echo "Converting bam to fastq ..." 
echo "**************************"

java -Xmx4g \
    -Djava.io.tmpdir=${work_dir}/tmp \
    -jar ${externaltoolsfolder}/SamToFastq.jar \
    INPUT=${work_dir}/${bam_name}.MT.bam \
    FASTQ=${work_dir}/${bam_name}.R1.fq \
    SECOND_END_FASTQ=${work_dir}/${bam_name}.R2.fq \
    UNPAIRED_FASTQ=${work_dir}/${bam_name}.fq \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=${work_dir}/tmp; 

echo "Done."; 
echo ""


echo "*****************************"
echo "Archiving processed bam files ..."
echo "*****************************"

mkdir -p ${work_dir}/processed_bam
mv ${work_dir}/*MT.bam ${work_dir}/*bai ${work_dir}/*sorted.bam ${work_dir}/processed_bam


echo "Done."
echo ""

echo "*************************"
echo "Aligning fastq to MT reference " ${reference}
echo "*************************"

bash align_mtDNA.sh ${reference} ${work_dir}/${bam_name} 

echo "Done."
echo ""


echo "*****************************"
echo "Removing unwanted  files ..."
echo "*****************************"
rm -rf ${work_dir}/processed_bam ${work_dir}/tmp ${work_dir}/*.fq ${work_dir}/*.unique_q10minimum.chrM.bam
echo "Done."
