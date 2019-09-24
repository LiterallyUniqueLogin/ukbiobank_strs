#!/bin/bash
#PBS -q hotel
#PBS -N convertToPfile
#PBS -l nodes=1:ppn=2 
#PBS -l walltime=00:10:00
#PBS -o /projects/ps-gymreklab/resources/ukbiobank/original/vcf_1_sample/hap/output/
#PBS -e /projects/ps-gymreklab/resources/ukbiobank/original/vcf_1_sample/hap/output/
#PBS -V
#PBS -M jmargoli@eng.ucsd.edu
#PBS -m a
#INPUT1 is the number of the chromosome to convert 

plink2 --pfile $UKB/original/pfile_converted/hap/chr${INPUT1} \
	--export vcf \
	--thin-indiv-count 1 \
	-out $UKB/original/vcf_1_sample/hap/chr$INPUT1
