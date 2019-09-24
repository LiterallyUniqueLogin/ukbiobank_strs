#!/bin/bash
#PBS -q hotel
#PBS -N convertToPfile
#PBS -l nodes=1:ppn=4 
#PBS -l walltime=00:20:00
#PBS -o /projects/ps-gymreklab/resources/ukbiobank/snpstr/pfile/conversion_output/
#PBS -e /projects/ps-gymreklab/resources/ukbiobank/snpstr/pfile/conversion_output/
#PBS -V
#PBS -M jmargoli@eng.ucsd.edu
#PBS -m a
##I only need two processors, but I need the memory that comes with four(?) processors
#INPUT1 is the number of the chromosome to convert 

plink2 --vcf $UKB/snpstr/1kg.snp.str.chr${INPUT1}.vcf.gz \
	--make-pgen \
	-out $UKB/snpstr/pfile/chr$INPUT1
