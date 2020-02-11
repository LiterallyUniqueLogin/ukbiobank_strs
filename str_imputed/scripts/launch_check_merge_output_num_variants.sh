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

sed -e "s/%RUN_NAME%/$1/g" \
	$UKB/str_imputed/check_merge_output_num_variants.pbs \
	> $TMPDIR/check_merge_output_num_variants_$1.pbs

source ~/.bashrc
conda activate bcftools
for start_pos in $(seq 1 5000000 250000000) ; do 
	end_pos=$((start_pos + 4999999))
	region_merge_file=$UKB/str_imputed/runs/$1/batches/chr${2}_pos_${start_pos}_to_${end_pos}.vcf.gz
	ls $region_merge_file
	reference_file=$UKB/snpstr/1kg.snp.str.chr$2.vcf.gz
	n_ref_variants=$(bcftools query -r $2:$start_pos-$end_pos -f '%POS\n' $reference_file | wc -l)
	echo "n_ref_variants for region $2:$start_pos-$end_pos is $n_ref_variants"
	qsub -v "INPUT1=$region_merge_file,INPUT2=$n_ref_variants" $TMPDIR/check_merge_output_num_variants_$1.pbs
done
conda deactivate
