source ~/.bashrc
conda activate bcftools
bcftools isec \
	$UKB/original/vcf_1_sample/hap/chr$1.vcf.gz \
        $UKB/snpstr/vcf_1_sample/chr$1.vcf.gz \
	-c none \
	-n=2 \
	#-o $UKB/pre_imputation_qc/common_variants/chr1.overlaps.txt

conda deactivate

