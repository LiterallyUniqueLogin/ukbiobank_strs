#!/bin/bash
#PBS -q hotel
#PBS -N convertToPfile
#PBS -l nodes=1:ppn=3 
#PBS -l walltime=00:10:00
#PBS -o /projects/ps-gymreklab/resources/ukbiobank/original/pfile_converted/hap/conversion_output/
#PBS -e /projects/ps-gymreklab/resources/ukbiobank/original/pfile_converted/hap/conversion_output/
#PBS -V
#PBS -M jmargoli@eng.ucsd.edu
#PBS -m a
##I only need two processors, but I need the memory that comes with three processors
#INPUT1 is the number of the chromosome to convert 

plink2 --bgen $UKB/original/bgen_original/hap/ukb_hap_chr${INPUT1}_v2.bgen ref-first \
	--sample $UKB/original/bgen_original/hap/*sample \
	--oxford-single-chr $INPUT1 \
	--make-pgen \
	-out $UKB/original/pfile_converted/hap/chr$INPUT1
