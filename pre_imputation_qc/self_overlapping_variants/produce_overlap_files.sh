#!/bin/bash

source ~/.bashrc
conda activate ukb_analysis
for chr in $(seq 1 22) ; do
	echo "Working on chr $chr"
	python $UKB/utilities/find_overlappying_variants.py \
		--vcf $UKB/snpstr/1kg.snp.str.chr${chr}.vcf.gz \
		--out $UKB/pre_imputation_qc/self_overlapping_variants/chr$chr \
		--filter $UKB/snpstr/duplicate_ids/chr*.txt
done
