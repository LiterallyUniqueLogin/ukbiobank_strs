#!/bin/bash

for i in {1..22}; do
	plink2 --pfile $UKB/snpstr/pfile/chr$i \
		--keep $UKB/snpstr/eur.sample \
		--freq \
		--out $UKB/pre_imputation_qc/allele_freqs/snpstr_eur/chr$i
done
