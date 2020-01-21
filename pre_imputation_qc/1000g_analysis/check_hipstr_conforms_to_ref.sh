#!/bin/bash

VCFS=""
for i in $(seq 1 22) ; do
	VCFS="$VCFS $DATASETS/1000Genomes/hipstr_calls_sample_subset/eur/chr$i/hipstr.chr$i.eur.vcf.gz"
done

source ~/.bashrc
conda activate ukb_analysis
python $UKB/utilities/check_vcf_against_ref.py \
	--ref $HUMAN/hg19/hg19.fa \
	--vcf $VCFS
conda deactivate
