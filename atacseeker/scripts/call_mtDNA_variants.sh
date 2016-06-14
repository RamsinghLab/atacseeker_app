## in_bam="chrM.CAP1727A2-041814-1-N.sorted.nodup"
in_bam=`dirname $1`"/"`basename $1 .bam`
externaltoolsfolder="./ext_tools/"
work_dir=`dirname ${in_bam}`
reference="../reference/chrRSRS.fa"
ref="RSRS"
UseIndelRealigner=true

echo "*****************************"
echo "Sorting, indexing and extraction of mitochondrial reads from bam file..." ${in_bam}.bam; 
echo "*****************************"
samtools sort ${in_bam}.bam ${in_bam}.sorted; 
samtools index ${in_bam}.sorted.bam; 
samtools view -b ${in_bam}.sorted.bam MT M chrMT chrM > ${in_bam}.MT.bam; 
echo "Done.";
echo ""

echo "*****************************"
echo "Converting bam to fastq..." ${in_bam}; #n=$(echo ${in_bam} | awk 'BEGIN{FS="."}{print $1}'); 
echo "*****************************"
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
echo ""

echo "*****************************"
echo "Archiving processed bam files"
echo "*****************************"
mkdir ${work_dir}/processed_bam
mv ${work_dir}/*MT.bam ${work_dir}/*bai ${work_dir}/*sorted.bam ${work_dir}/processed_bam
tar -czf ${work_dir}/processed_bam.tar.gz ${work_dir}/processed_bam
#rm -r ${work_dir}/processed_bam

echo "*************************"
echo "Aligning fastq to chrRSRS"
echo "*************************"
bash align_mtDNA.sh ${reference} ${in_bam}
echo "Done."
echo ""

echo "*******************************"
echo "archiving processed fastq files"
echo "*******************************"
mkdir ${work_dir}/processed_fastq
mv ${in_bam}*fastq ${work_dir}/processed_fastq
tar czf ${work_dir}/processed_fastq.tar.gz ${work_dir}/processed_fastq
#rm -r ${work_dir}/processed_fastq 


GENOME_TYPE="chrRSRS" 
OUT_BAM=${in_bam}".unique_q10minimum."$GENOME_TYPE".nodup.bam"
SORTED=`basename $OUT_BAM .bam`.sorted
DIR=`dirname $OUT_BAM`

if $UseIndelRealigner
	then
		echo "*********************************************************"
		echo "* REALIGNING KNOWN INDELS WITH GATK INDELREALIGNER..."
		echo "* Realigning known indels for file" ${SORTED}".bam using data/MITOMAP_HMTDB_known_indels.vcf as reference..."
		echo "*********************************************************"
		java -Xmx4g \
			-Djava.io.tmpdir=`pwd`/tmp \
			-jar ${externaltoolsfolder}GenomeAnalysisTK.jar \
			-T IndelRealigner \
			-R ./data/chr${ref}.fa \
			-I $DIR"/"${SORTED}".bam" \
			-o $DIR"/"${SORTED}".realigned.bam" \
			-targetIntervals ./data/intervals_file_${ref}.list  \
			-known ./data/MITOMAP_HMTDB_known_indels_${ref}.vcf \
			-compress 0;
		echo "Done."
		echo ""
		SORTED=${SORTED}".realigned"
fi

echo "**********************************************"
echo "** removing duplicates and output SAM file"
echo "**********************************************"
samtools rmdup $DIR"/"${SORTED}".bam" $DIR"/"${SORTED}".nodup.bam"
samtools index $DIR"/"${SORTED}".nodup.bam"
samtools view $DIR"/"${SORTED}".nodup.bam" > $DIR"/"${SORTED}".nodup.sam"
echo "Done."
echo ""

echo "********************************"
echo "*** assemble mtDNA"
echo "********************************"
python assembleMTgenome.py \
	-r "../reference/" \
	-f chrRSRS.fa \
	-i $DIR"/"${SORTED}".nodup.bam" \
	-a hg19.fa \
	-s `which samtools` -FCUP \
	-o $DIR"/"${in_bam}
echo "Done."
echo ""

echo "********************************"
echo "***  Call mtDNA Variants"
echo "********************************"
python VCFoutput.py -r RSRS -i $DIR
echo "Done."
echo ""