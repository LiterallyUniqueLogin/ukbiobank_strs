#!/bin/bash

for i in {1..22}; do
	plink2 --pfile $UKB/original/pfile_converted/hap/chr$i \
		--keep $UKB/sample_qc/snp_informed/hap/combined/white_brit_high_quality_sane_kinship_unrelated.sample \
		--freq \
		--out $UKB/pre_imputation_qc/allele_freqs/hap_white_brit_high_quality_sane_kinship_unrelated/chr$i
done
