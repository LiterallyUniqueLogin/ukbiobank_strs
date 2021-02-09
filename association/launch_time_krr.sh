#!/usr/bin/env bash
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=7:00:00
#SBATCH --job-name=time_krr
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbmargoliash@ucsd.edu
#SBATCH --output=/expanse/projects/gymreklab/jmargoli/ukbiobank/association/output/%x_%j.out
#SBATCH --account=sds154

source "$HOME"/.bash_profile
conda activate ukb

cd "$UKB/association" || { echo "Directory $UKB/association doesn't exist" ; exit 1 ; }

python time_krr.py

conda deactivate
