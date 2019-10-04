#!/bin/bash

source ~/.bashrc

#Inputs one and two are the beginning and end numbers of the samples to process
#counting starting at 1, inclusive.
#Input 3 is the chrom number
#E.g. 2001-4000 being the second group of 2000 samples
startLine=$(($1+2))
endLine=$(($2+2))
exitLine=$(($2+3))

TMP_SAMPLE_FILE=$TMPDIR/samples_$1_to_$2.sample
TMP_INPUT_VCF_NOEXT=$TMPDIR/chr$3_samples_$1_to_$2
sed -n "1p;${startLine},${endLine}p;${exitLine}q" \
 $UKB/original/bgen_original/hap/*sample \
 > $TMP_SAMPLE_FILE

plink2 --pfile $UKB/original/pfile_converted/hap/chr$3 \
	--export vcf \
	--out $TMP_INPUT_VCF_NOEXT \
	--ref-allele force $UKB/non_genetic_data/showcase/ukb_snp_bim/ukb_snp_chr$3_v2.bim 5 2 \
	--keep $TMP_SAMPLE_FILE

conda activate bcftools
bgzip $TMP_INPUT_VCF_NOEXT.vcf

tabix $TMP_INPUT_VCF_NOEXT.vcf.gz
conda deactivate
