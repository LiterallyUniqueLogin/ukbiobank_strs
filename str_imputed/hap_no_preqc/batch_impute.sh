#!/bin/bash

#first argument is the chromosome number to do the imputation for
if [ -z "$1" ] ; then 
    echo "Error: Need first argument - the chromosome number"
    exit -1
fi

numSamples=487409
for (( minId=1 ; minId<=numSamples ; minId+=1000 )) ; do
	maxId=$((minId+999))
	maxId=$(( maxId<numSamples ? maxId : numSamples )) #Take the smaller of the two bounds

	OUTFILE=$UKB/str_imputed/hap_no_preqc/vcf_batches/chr${1}_samples_${minId}_to_${maxId}.vcf.gz
	if [ -f "$OUTFILE" ] ; then 
		echo $OUTFILE already exists. Skipping this batch.
		continue
	fi

	qsub -v "INPUT1=$minId,INPUT2=$maxId,INPUT3=$1" batch_script.pbs
done


