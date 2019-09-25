#!/bin/bash
#PBS -q hotel
#PBS -N snpstr_snpstr_intersect
#PBS -l nodes=1:ppn=1
#PBS -l walltime=03:00:00
#PBS -o /projects/ps-gymreklab/resources/ukbiobank/pre_imputation_qc/overlapping_variants/snpstr_snpstr/output/
#PBS -e /projects/ps-gymreklab/resources/ukbiobank/pre_imputation_qc/overlapping_variants/snpstr_snpstr/output/
#PBS -V
#PBS -M jmargoli@eng.ucsd.edu
#PBS -m a
#INPUT1 is the number of the chromosome to convert 

source ~/.bashrc
conda activate bedtools
echo "Running bedtools"
time bedtools intersect -wa -wb \
	-a $UKB/snpstr/vcf_1_sample/chr${INPUT1}.vcf \
	-b $UKB/snpstr/vcf_1_sample/chr${INPUT1}.vcf \
	-sorted \
	> $TMPDIR/all_overlap.txt
echo "Cutting first overlap file"
time cut -f1,2,4,5 $TMPDIR/all_overlap.txt > $TMPDIR/overlap_1.txt
echo "Cutting second overlap file"
time cut -f11,12,14,15 $TMPDIR/all_overlap.txt > $TMPDIR/overlap_2.txt
echo "Diffing overlap files"
time diff $TMPDIR/overlap_1.txt $TMPDIR/overlap_2.txt \
	> $UKB/pre_imputation_qc/overlapping_variants/snpstr_snpstr/chr${INPUT1}.txt

