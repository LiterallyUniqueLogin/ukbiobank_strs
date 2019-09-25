#!/bin/bash
#PBS -q hotel
#PBS -N beagle_test
#PBS -l nodes=1:ppn=14
#PBS -l walltime=01:00:00
#PBS -o /projects/ps-gymreklab/resources/ukbiobank/pre_imputation_qc/beagle_check/std.out
#PBS -e /projects/ps-gymreklab/resources/ukbiobank/pre_imputation_qc/beagle_check/err.out
#PBS -V
#PBS -M jmargoli@eng.ucsd.edu
#PBS -m ae
#PBS -l pmem=80gb
#Take as many cores as possible (tscc max is 14/node) (if I recall, Beagle used ~16 cores on Snorlax)

source ~/.bashrc
conda activate java8
java -Xmx50000m -jar $SOURCE/beagle.24Aug19.3e8.jar \
	ref=$UKB/snpstr/1kg.snp.str.chr9.vcf.gz \
	gt=$UKB/pre_imputation_qc/beagle_check/chr9_samples_1_to_2000_outof_487409.vcf \
	out=$UKB/pre_imputation_qc/beagle_check/imputed \
	map=$UKB/genetic_maps/plink.chr9.GRCh37.map \
	impute=true gp=true ap=true
