#!/bin/bash

thresh=0.0001
for i in {1..22}; do
        plink2 --pfile $UKB/original/pfile_converted/hap/chr$i \
		--maf $thresh \
		--write-snplist \
                --out $UKB/pre_imputation_qc/allele_freqs/hap/maf_filter/min_maf_$thresh/chr$i
done
