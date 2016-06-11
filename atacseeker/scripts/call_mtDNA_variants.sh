in_bam="chrM.CAP1727A2-041814-1-N.sorted.nodup"
externaltoolsfolder="/Users/asifzubair/projects/MToolBox/MToolBox/atac_3vs3_data/ext_tools/"
work_dir=`pwd`
reference="../reference/chrRSRS.fa"
ref="RSRS"
UseIndelRealigner=false


echo "Sorting, indexing and extraction of mitochondrial reads from bam file..." ${in_bam}.bam; 
samtools sort ${in_bam}.bam ${in_bam}.sorted; 
samtools index ${in_bam}.sorted.bam; 
samtools view -b ${in_bam}.sorted.bam MT M chrMT chrM > ${in_bam}.MT.bam; 
echo "Done.";
echo ""


echo "Converting bam to fastq..." ${in_bam}; #n=$(echo ${in_bam} | awk 'BEGIN{FS="."}{print $1}'); 
java -Xmx4g \
	-Djava.io.tmpdir=${work_dir}/tmp \
	-jar ${externaltoolsfolder}SamToFastq.jar \
	INPUT=${in_bam}.MT.bam \
	FASTQ=${in_bam}.R1.fastq \
	SECOND_END_FASTQ=${in_bam}.R2.fastq \
	UNPAIRED_FASTQ=${in_bam}.fastq \
	VALIDATION_STRINGENCY=SILENT \
	TMP_DIR=${work_dir}/tmp; 
echo "Done."; 


echo "Archiving processed bam files"
mkdir ${work_dir}/processed_bam
mv ${work_dir}/*MT.bam ${work_dir}/*bai ${work_dir}/*sorted.bam ${work_dir}/processed_bam
tar -czf ${work_dir}/processed_bam.tar.gz ${work_dir}/processed_bam
#rm -r ${work_dir}/processed_bam


echo "Aligning fastq to chrRSRS"
bash align_mtDNA.sh ${reference} ${in_bam}


echo "archiving processed fastq files"
mkdir ${work_dir}/processed_fastq
mv ${in_bam}*fastq ${work_dir}/processed_fastq
tar czf ${work_dir}/processed_fastq.tar.gz ${work_dir}/processed_fastq
#rm -r ${work_dir}/processed_fastq 


GENOME_TYPE="chrRSRS" 
OUT_BAM=${in_bam}".unique_q10minimum."$GENOME_TYPE".nodup.bam"
SORTED=`basename $OUT_BAM .bam`.sorted


if $UseIndelRealigner
	then
		echo ""
		echo "##### REALIGNING KNOWN INDELS WITH GATK INDELREALIGNER..."
		echo ""
		echo "Realigning known indels for file" ${SORTED}".bam using data/MITOMAP_HMTDB_known_indels.vcf as reference..."
		java -Xmx4g \
			-Djava.io.tmpdir=`pwd`/tmp \
			-jar ${externaltoolsfolder}GenomeAnalysisTK.jar \
			-T IndelRealigner \
			-R ./data/chr${ref}.fa \
			-I ${SORTED}".bam" \
			-o ${SORTED}".realigned.bam" \
			-targetIntervals ./data/intervals_file_${ref}.list  \
			-known ./data/MITOMAP_HMTDB_known_indels_${ref}.vcf \
			-compress 0;
		SORTED=${SORTED}".realigned"
fi


echo "removing duplicates and output SAM file"
samtools rmdup ${SORTED}".bam" ${SORTED}".nodup.bam"
samtools index ${SORTED}".nodup.bam"
samtools view ${SORTED}".nodup.bam" > ${SORTED}".nodup.sam"


echo "assemble mtDNA"
assembleMTgenome.py \
	-r "../reference/" \
	-f chrRSRS.fa \
	-i ${SORTED}".nodup.bam" \
	-a hg19.fa \
	-s `which samtools` -FCUP \
	-o ${in_bam}
