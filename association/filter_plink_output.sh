#!/usr/bin/env bash
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=7:30:00
#SBATCH --job-name=filter_plink_output
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbmargoliash@ucsd.edu
#SBATCH --output=/expanse/projects/gymreklab/jmargoli/ukbiobank/association/output/%x_%j.out
#SBATCH --account=sds154

dir=/expanse/projects/gymreklab/jmargoli/ukbiobank/association/runs/dosages_snps_21_plink/results/chr21
for phenotype in height total_bilirubin ; do
	grep ADD $dir/plink2.${phenotype}_inv_norm_rank.glm.linear \
		> $dir/plink2.${phenotype}_inv_norm_rank.glm.linear.summary
done
