#!/bin/bash

for i in {1..22}; do
	plink2 --pfile $UKB/original/pfile_converted/hap/chr$i \
		--freq \
		--out $UKB/pre_imputation_qc/allele_freqs/hap/chr$i
done
