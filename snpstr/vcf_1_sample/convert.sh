#!/bin/bash
#PBS -q hotel
#PBS -N convertToPfile
#PBS -l nodes=1:ppn=2 
#PBS -l walltime=00:10:00
#PBS -o /projects/ps-gymreklab/resources/ukbiobank/snpstr/vcf_1_sample/output/
#PBS -e /projects/ps-gymreklab/resources/ukbiobank/snpstr/vcf_1_sample/output/
#PBS -V
#PBS -M jmargoli@eng.ucsd.edu
#PBS -m a
#INPUT1 is the number of the chromosome to convert 

plink2 --pfile $UKB/snpstr/pfile/chr${INPUT1} \
	--export vcf ref-first \
	--thin-indiv-count 1 \
	-out $UKB/snpstr/vcf_1_sample/chr$INPUT1
	


