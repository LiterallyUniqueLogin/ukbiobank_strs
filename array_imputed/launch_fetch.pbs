#!/bin/bash
#PBS -q hotel
#PBS -N fetch_imputed_files 
#PBS -l nodes=1:ppn=1
#PBS -l walltime=96:00:00
#PBS -o /projects/ps-gymreklab/jmargoli/ukbiobank/array_imputed/output
#PBS -e /projects/ps-gymreklab/jmargoli/ukbiobank/array_imputed/output
#PBS -V
#PBS -M jmargoli@ucsd.edu
#PBS -m a

source "$HOME"/.bashrc
cd "$UKB/array_imputed" || { echo "Directory $UKB/array_imputed doesn't exist" ; exit 1 ; }
conda activate ukb
python fetch_imputed_files.py 
conda deactivate

