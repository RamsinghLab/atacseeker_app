for s in 071114-3 100214-1 102414-3 103114-1 103114-2; do
	bash call_mtDNA_variants.sh /data/chrM.${s}_N.sorted.nodup.bam RSRS
	bash call_mtDNA_variants.sh /data/chrM.${s}_S.sorted.nodup.bam RSRS
	python VCFoutput.py RSRS
	mv VCF_file.vcf ${s}_NvsS.vcf
	mv VCF_dict_tmp VCF_dict_${s}
done
