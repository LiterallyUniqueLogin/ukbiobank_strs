#!/bin/bash

source ~/.bashrc

#Inputs one and two are the beginning and end numbers of the samples to process
#counting starting at 1, inclusive.
#E.g. 2001-4000 being the second group of 2000 samples
#Input 3 is the chrom number
#Input 4 is the directory containing the pfiles
#Input 5 is the .sample file
#This file must have the format 
#(or any other format that plink --keep will recognize)
#ID_1 ID_2 missing sex
#2497795 2497795 0 1
#2143467 2143467 0 2
#... 
#Input 6 must be the run name

startLine=$(($1+2))
endLine=$(($2+2))
exitLine=$(($2+3))

mkdir $TMPDIR/$6
TMP_SAMPLE_FILE=$TMPDIR/$6/samples_$1_to_$2.sample
TMP_INPUT_VCF_NOEXT=$TMPDIR/$6/chr$3_samples_$1_to_$2
sed -n "1p;${startLine},${endLine}p;${exitLine}q" $5/*.sample > $TMP_SAMPLE_FILE

plink2 --pfile $4/chr$3 \
	--export vcf \
	--out $TMP_INPUT_VCF_NOEXT \
	--keep $TMP_SAMPLE_FILE

conda activate bcftools
bgzip $TMP_INPUT_VCF_NOEXT.vcf

tabix $TMP_INPUT_VCF_NOEXT.vcf.gz
conda deactivate
