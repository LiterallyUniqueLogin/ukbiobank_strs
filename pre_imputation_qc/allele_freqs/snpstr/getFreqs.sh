#!/bin/bash

for i in {1..22}; do
	plink2 --pfile $UKB/snpstr/pfile/chr$i \
		--freq \
		--out $UKB/pre_imputation_qc/allele_freqs/snpstr/chr$i
done
