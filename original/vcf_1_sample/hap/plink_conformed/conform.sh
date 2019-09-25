#!/bin/bash

INPUT1=22
plink2 --vcf $UKB/original/vcf_1_sample/hap/chr${INPUT1}.vcf \
	--ref-allele force $UKB/snpstr/vcf_1_sample/chr${INPUT1}.vcf 4 3 '#' \
	--export vcf \
	--out $UKB/original/vcf_1_sample/hap/plink_conformed/chr${INPUT1}
