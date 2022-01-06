#!/usr/bin/env bash

outfname=$1
ethnicity=$2
phenotype=$3
chrom=$4
start_pos=$5
end_pos=$6
str_pos=$7
load_vars=$8

if [[ ("$#" -lt 7) || ("$#" -gt 8) ]] ; then
	echo "Illegal number of arguments. Expected 6 or 7 (outfname, ethnicity, phenotype, chrom, start_pos, end_pos, str_pos, load_vars)." 1>&2
	exit 1
fi

#plink2_prefix="$UKB"/microarray/pfile_converted/chr${chrom}
plink2_prefix="$UKB"/array_imputed/pfile_converted/chr${chrom}
plink2_pgen=${plink2_prefix}.pgen
if [[ ! -f "${plink2_prefix}.pgen" ]] ; then
	echo "pgen file ${plink2_pgen} doesn't exist" 1>&2
	exit 2
fi

sample_file="$UKB"/sample_qc/runs/${ethnicity}/${phenotype}/combined_unrelated.sample
if [[ ! -f ${sample_file} ]] ; then
	echo "sample file ${sample_file} doesn't exit" 1>&2
	exit 3
fi

run_name="${ethnicity}_${phenotype}_${chrom}_${start_pos}_${end_pos}_${str_pos}"
plink1_dir="post_finemapping/intermediate_results/plink1_chrom_regions/${run_name}"
mkdir -p "$plink1_dir"
cd "$plink1_dir" || {
	echo 'cd failed'
	exit 5	
}

sample_tmpfile=$(mktemp --tmpdir XXXXXXXXXX.sample)
{ echo '#FID IID' ; tail -n +2 "${sample_file}" | awk '{print $1 " " $1 }' ; } > "$sample_tmpfile"
plink2_cmd="
$UKB/utilities/plink2 \
       	--pfile \"${plink2_prefix}\"
	--keep \"$sample_tmpfile\"
	--snps-only
	--chr \"$chrom\" --from-bp \"$start_pos\" --to-bp \"$end_pos\"
	--make-bed
	--seed 13 --thin-indiv-count 10000
"
if [[ "$load_vars" == "" ]] ; then
	plink2_cmd="$plink2_cmd --mac 1"
fi
echo "plink2_cmd: $plink2_cmd"
eval $plink2_cmd

if [[ "$load_vars" != "" ]] ; then
	if [[ "$(wc -w "$load_vars")" -ne "$(wc -l plink2.bim)" ]] ; then
		echo "Did not load correct number of variants" 2>&1
		exit  6
	fi
fi

rm "$sample_tmpfile"
rm plink2.log

for suffix in bed bim fam ; do
	mv plink2.${suffix} "data.${suffix}"
done

$UKB/utilities/plink --bfile data --r2 square0 gz 
mv plink.ld.gz ld.gz
rm plink.log

cd $UKB
cmd="$UKB/post_finemapping/my_haploview.py $plink1_dir/ld.gz $plink1_dir/data.bim $outfname $str_pos"
if [[ "$load_vars" == "" ]] ; then
	cmd="$cmd --savevars $plink1_dir/var_subset.txt"
else
	cmd="$cmd --loadvars post_finemapping/intermediate_results/plink1_chrom_regions/white_brits_${phenotype}_${chrom}_${start_pos}_${end_pos}_${str_pos}/var_subset.txt"
fi

eval $cmd


