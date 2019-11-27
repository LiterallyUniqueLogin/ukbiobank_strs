#!/bin/bash
for chr in $(seq 1 22) ; do 
	mkdir $UKB/pre_imputation_qc/ref_panel_size/output/variant_overlaps/chr$chr

	bcftools isec -n+2 -c all \
		-p $UKB/pre_imputation_qc/ref_panel_size/output/variant_overlaps/chr$chr \
		$UKB/snpstr/vcf_1_sample/chr$chr.vcf.gz \
		$DATASETS/1000Genomes/hipstr_calls_sample_subset/eur/chr$chr/hipstr.chr$chr.eur.vcf.gz
done
