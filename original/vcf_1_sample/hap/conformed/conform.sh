#!/bin/bash

INPUT1=22
java -Xmx8g -jar $SOURCE/conform-gt.24May16.cee.jar \
	ref=$UKB/snpstr/vcf_1_sample/chr${INPUT1}.vcf \
	gt=$UKB/original/vcf_1_sample/hap/chr${INPUT1}.vcf \
	chrom=$INPUT1 \
	match=POS \
	out=$UKB/original/vcf_1_sample/hap/conformed/chr${INPUT1}
