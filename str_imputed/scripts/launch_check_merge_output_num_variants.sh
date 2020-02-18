#!/bin/bash

if [ -z "$1" ] ; then
        echo "didn't give first argument - should be the run name" 1>&2
        exit -1
fi

if [ -z "$2" ] ; then
        echo "didn't give second argument - should be chrom number" 1>&2
        exit -1
fi

#if more arguments are supplied
#then they should be regions and only those regions will be launched

if [ -z "$TMPDIR" ] ; then
        echo "Didn't set TMPDIR" 1>&2
        exit -1
fi


source ~/.bashrc
conda activate bcftools
for start_pos in $(seq 1 5000000 250000000) ; do 
	if (($# >= 3)) ; then
		SKIP="yep"
		for (( idx=3 ; idx<=$# ; idx++)) ; do 
			if [ "${!idx}" == "$start_pos" ] ; then
				SKIP=""
				break
			fi
		done
		if [ "$SKIP" == "yep" ] ; then
			continue
		fi
	fi

	sed -e "s/%RUN_NAME%/$1/g" -e "s/%CHROM%/$2/g" -e "s/%POS%/$start_pos/g"\
		$UKB/str_imputed/scripts/check_merge_output_num_variants.pbs \
		> $TMPDIR/check_merge_output_num_variants_${1}_${2}_${start_pos}.pbs

	end_pos=$((start_pos + 4999999))
	region_merge_file=$UKB/str_imputed/runs/$1/batches/chr${2}_pos_${start_pos}_to_${end_pos}.vcf.gz
	reference_file=$UKB/snpstr/1kg.snp.str.chr$2.vcf.gz
	n_ref_variants=$(bcftools query -r $2:$start_pos-$end_pos -f '%POS\n' $reference_file | wc -l)
	echo "n_ref_variants for region $2:$start_pos-$end_pos is $n_ref_variants"
	qsub -v "INPUT1=$region_merge_file,INPUT2=$n_ref_variants" $TMPDIR/check_merge_output_num_variants_${1}_${2}_${start_pos}.pbs
done
conda deactivate
