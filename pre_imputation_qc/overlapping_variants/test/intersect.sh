TMPDIR=./temp
mkdir $TMPDIR

source ~/.bashrc
conda activate bedtools
echo "Running bedtools"
time bedtools intersect -wa -wb \
	-a $UKB/pre_imputation_qc/overlapping_variants/test/test1.vcf \
	-b $UKB/pre_imputation_qc/overlapping_variants/test/test2.vcf \
	> $TMPDIR/all_overlap.txt
echo "Cutting first overlap file"
time cut -f1,2,4,5 $TMPDIR/all_overlap.txt > $TMPDIR/overlap_1.txt
echo "Cutting second overlap file"
time cut -f11,12,14,15 $TMPDIR/all_overlap.txt > $TMPDIR/overlap_2.txt
echo "Diffing overlap files"
time diff $TMPDIR/overlap_1.txt $TMPDIR/overlap_2.txt \
	> $UKB/pre_imputation_qc/overlapping_variants/test/test_overlap.txt

conda deactivate
