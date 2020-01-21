#!/bin/bash

source ~/.bashrc
conda activate ukb_analysis
python $UKB/utilities/check_vcf_against_ref.py \
	--ref $HUMAN/hg19/hg19.fa \
	--vcf $UKB/snpstr/1kg.snp.str.chr*.vcf.gz
conda deactivate
