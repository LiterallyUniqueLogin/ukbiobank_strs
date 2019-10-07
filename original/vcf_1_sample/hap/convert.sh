plink2 --pfile $UKB/original/pfile_converted/hap/chr${INPUT1} \
	--export vcf ref-first \
	--thin-indiv-count 1 \
	--ref-allele force $UKB/non_genetic_data/showcase/ukb_snp_bim/ukb_snp_chr${INPUT1}_v2.bim 5 2 \
	-out $UKB/original/vcf_1_sample/hap/chr$INPUT1

conda activate	

