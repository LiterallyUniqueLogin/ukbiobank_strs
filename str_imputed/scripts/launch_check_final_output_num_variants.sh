#!/bin/bash

if [ -z "$1" ] ; then
        echo "didn't give first argument - should be the run name" 1>&2
        exit -1
fi

if [ -z "$2" ] ; then
        echo "didn't give second argument - should be chrom number" 1>&2
        exit -1
fi

if [ -z "$TMPDIR" ] ; then
        echo "Didn't set TMPDIR" 1>&2
        exit -1
fi


sed -e "s/%RUN_NAME%/$1/g" -e "s/%CHROM%/$2/g" \
	$UKB/str_imputed/scripts/check_final_output_num_variants.pbs \
	> $TMPDIR/check_final_output_num_variants_${1}_${2}.pbs

source ~/.bashrc
conda activate bcftools

reference_file=$UKB/snpstr/1kg.snp.str.chr$2.vcf.gz
n_ref_variants=$(bcftools query -f '%POS\n' $reference_file | wc -l)

conda deactivate

echo "n_ref_variants for chromosome $2 is $n_ref_variants"
qsub -v "INPUT1=$n_ref_variants" $TMPDIR/check_final_output_num_variants_${1}_${2}.pbs

