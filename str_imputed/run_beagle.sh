#!/bin/bash

source ~/.bashrc

if [ -z "$1" ] ; then
        echo "didn't give first argument - should be the first sample number to process \
(1 indexed, inclusive)"
        exit -1
fi

if [ -z "$2" ] ; then
        echo "didn't give second argument - should be the last sample number to process \
(1 indexed, inclusive)"
        exit -1
fi

if [ -z "$3" ] ; then
        echo "didn't give third argument - should be the chromosome number"
        exit -1
if

if [ -z "$4" ] ; then
        echo "didn't give fourth argument - should be the run name"
        exit -1
fi

if [ -z "$TMPDIR" ] ; then
        echo "Didn't set TMPDIR"
        exit -1
fi

TMP_INPUT_VCF_NOEXT=$TMPDIR/$4/chr$3_samples_$1_to_$2
OUTFILE_NOEXT=$UKB/str_imputed/runs/$4/batches/chr$3_samples_$1_to_$2 

conda activate java8
java -Xmx12500m -jar $SOURCE/beagle.25Nov19.28d.jar \
	ref=$UKB/snpstr/1kg.snp.str.chr$3.vcf.gz \
	gt=$TMP_INPUT_VCF_NOEXT.vcf.gz \
	out=$OUTFILE_NOEXT \
	map=$SOURCE/beagle_genetic_maps/plink.chr$3.GRCh37.map \
	impute=true gp=false ap=true
conda deactivate

conda activate bcftools
tabix $OUTFILE_NOEXT.vcf.gz
conda deactivate

