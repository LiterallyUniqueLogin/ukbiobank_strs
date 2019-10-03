source ~/.bashrc
conda activate bcftools
bgzip -c $UKB/original/vcf_1_sample/hap/chr${INPUT1}.vcf > \
	 $UKB/original/vcf_1_sample/hap/chr${INPUT1}.vcf.gz
conda deactivate
conda activate bedtools
echo "Running bedtools"
time bedtools intersect -wa -wb \
	-a $UKB/snpstr/vcf_1_sample/chr${INPUT1}.vcf \
	-b $UKB/original/vcf_1_sample/hap/chr${INPUT1}.vcf.gz \
	-sorted \
	> $TMPDIR/all_overlap.txt
echo "Cutting first overlap file"
time cut -f1,2,4,5 $TMPDIR/all_overlap.txt > $TMPDIR/overlap_1.txt
echo "Cutting second overlap file"
time cut -f11,12,14,15 $TMPDIR/all_overlap.txt > $TMPDIR/overlap_2.txt
echo "Diffing overlap files"
time diff $TMPDIR/overlap_1.txt $TMPDIR/overlap_2.txt \
	> $UKB/pre_imputation_qc/overlapping_variants/snpstr_hap/chr${INPUT1}.txt
conda deactivate
