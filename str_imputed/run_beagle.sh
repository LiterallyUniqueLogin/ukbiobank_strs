#!/bin/bash

source ~/.bashrc

#Inputs one and two are the beginning and end numbers of the samples to process
#counting starting at 1, inclusive.
#Input 3 is the chrom number
#E.g. 2001-4000 being the second group of 2000 samples
#Input 4 is the name of the run dir

TMP_INPUT_VCF_NOEXT=$TMPDIR/$4/chr$3_samples_$1_to_$2
OUTFILE_NOEXT=$UKB/str_imputed/$4/batches/chr$3_samples_$1_to_$2 

conda activate java8
java -Xmx12500m -jar $SOURCE/beagle.21Sep19.ec3.jar \
	ref=$UKB/snpstr/1kg.snp.str.chr$3.vcf.gz \
	gt=$TMP_INPUT_VCF_NOEXT.vcf.gz \
	out=$OUTFILE_NOEXT \
	map=$SOURCE/beagle_genetic_maps/plink.chr$3.GRCh37.map \
	impute=true gp=false ap=true
conda deactivate

#create both the .vcf file and .vcf.gz file because
#I haven't finalized how I want to do merging
conda activate bcftools
bgzip -d $OUTFILE_NOEXT.vcf.gz
$UKB/utilities/add_ap_metadata.sh $OUTFILE_NOEXT no-gz
bgzip -c $OUTFILE_NOEXT.vcf > $OUTFILE_NOEXT.vcf.gz
tabix $OUTFILE_NOEXT.vcf.gz
conda deactivate
