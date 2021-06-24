#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --account=ddp268
#SBATCH --job-name=launch_snakemake
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=47:30:00
#SBATCH --output=/expanse/projects/gymreklab/jmargoli/ukbiobank/workflow/logs/%x_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jonathan.margoliash@gmail.com
#SBATCH --export=ALL

if [ -z "$TARGET" ] ; then echo "TARGET variable is unset" ; exit 1 ; fi
echo "TARGET $TARGET"
echo "TARGET $TARGET" >&2

# I've installed snakemake with conda, so I need to start
# anaconda before calling snakemake
source "$HOME"/.bashrc
conda activate ukb
cd "$UKB" || { echo "Directory $UKB doesn't exist" ; exit 1 ; }

# intentional no quote on TARGET
snakemake \
	--profile workflow \
	$TARGET

conda deactivate

