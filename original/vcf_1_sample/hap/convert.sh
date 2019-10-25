plink2 --pfile $UKB/original/pfile_converted/hap/chr${INPUT1} \
	--export vcf ref-first \
	--thin-indiv-count 1 \
	-out $UKB/original/vcf_1_sample/hap/chr$INPUT1

conda activate	

