#!/bin/bash
#PBS -q hotel
#PBS -N snpstr_snpstr isec
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:15:00
#PBS -o /projects/ps-gymreklab/resources/ukbiobank/pre_imputation_qc/overlapping_variants/snpstr_snpstr/output/
#PBS -e /projects/ps-gymreklab/resources/ukbiobank/pre_imputation_qc/overlapping_variants/snpstr_snpstr/output/
#PBS -V
#PBS -M jmargoli@eng.ucsd.edu
#PBS -m a
#INPUT1 is the number of the chromosome to convert 

source ~/.bashrc
conda activate bedtools
bedtools intersect -wa -wb \
	-a $UKB/snpstr/1kg.snp.str.chr${INPUT1}.vcf.gz \
	-b $UKB/snpstr/1kg.snp.str.chr${INPUT1}.vcf.gz \
	-sorted \
	> $TMPDIR/all_overlap.txt
cut -f1,2,4,5 $TMPDIR/all_overlap.txt > $TMPDIR/overlap_1.txt
cut -f10,11,13,14 $TMPDIR/all_overlap.txt > $TMPDIR/overlap_2.txt
diff $TMPDIR/overlap_1.txt $TMPDIR/overlap_2.txt \
	> $UKB/pre_imputation_qc/overlapping_variants/snpstr_snpstr/chr${INPUT1}.txt

