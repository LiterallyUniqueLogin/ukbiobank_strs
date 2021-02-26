#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --account=sds154
#SBATCH --job-name=launch_linear_association
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=5G
#SBATCH --time=47:30:00
#SBATCH --output=/expanse/projects/gymreklab/jmargoli/ukbiobank/association/output/%x_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbmargoliash@ucsd.edu
#SBATCH --export=ALL

if [ -z "$ARGS" ] ; then echo "ARGS variable is unset" ; exit 1 ; fi
echo "ARGS $ARGS"
echo "ARGS $ARGS" >&2

source "$HOME"/.bashrc
cd "$UKB/association" || { echo "Directory $UKB/association doesn't exist" ; exit 1 ; }
conda activate ukb

# intentional no quote here
python linear_association.py $ARGS

conda deactivate

