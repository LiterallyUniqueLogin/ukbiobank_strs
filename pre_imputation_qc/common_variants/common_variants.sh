dir=$UKB/pre_imputation_qc/common_variants
# use pfiles since they have been ref/alt recoded
cut -f1,2,4,5 -d $'\t' $UKB/microarray/pfile_converted/chr*.pvar > $dir/hap_vars.tab
{ for chrom in $(seq 1 22) ; do bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" $UKB/snpstr/1kg.snp.str.chr${chrom}.vcf.gz ; done } > $dir/snpstr_vars.tab
comm -12 <(sed 's/[a-z]/\U&/g' $dir/hap_vars.tab | sort) <(sed 's/[a-z]/\U&/g' $dir/snpstr_vars.tab | sort) > $dir/common_vars.tab
