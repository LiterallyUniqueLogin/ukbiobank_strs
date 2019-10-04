#!/bin/bash

source ~/.bashrc

#Inputs one and two are the beginning and end numbers of the samples to process
#counting starting at 1, inclusive.
#Input 3 is the chrom number
#E.g. 2001-4000 being the second group of 2000 samples

TMP_INPUT_VCF_NOEXT=$TMPDIR/chr$3_samples_$1_to_$2

conda activate java8
java -Xmx10000m -jar $SOURCE/beagle.21Sep19.ec3.jar \
	ref=$UKB/snpstr/1kg.snp.str.chr$3.vcf.gz \
	gt=$TMP_INPUT_VCF_NOEXT.vcf.gz \
	out=$UKB/str_imputed/hap_no_preqc/vcf_batches/chr$3_samples_$1_to_$2 \
	map=$SOURCE/beagle_genetic_maps/plink.chr$3.GRCh37.map \
	impute=true gp=false ap=true
conda deactivate
