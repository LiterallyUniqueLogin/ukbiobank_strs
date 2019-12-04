#!/bin/bash

EUR_HIPSTR_AND_REF_SAMPLES=$UKB/pre_imputation_qc/ref_panel_size/output/compare_samples/eur_hipstr_and_ref.sample
for sample in $(tail -n+2 $EUR_HIPSTR_AND_REF_SAMPLES) ; do
	echo $sample
	qsub -v "INPUT1=$sample" compare_imputations.pbs
done

