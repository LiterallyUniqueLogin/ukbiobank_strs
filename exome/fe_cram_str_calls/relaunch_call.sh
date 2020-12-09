#!/usr/bin/env bash

if [ -z "$1" ] ; then
	echo "Didn't give input STR" 
	exit 1
fi
if [ ! -f "str_beds/$1.bed" ] ; then
	echo "str_beds/$1.bed doesn't exit"
	exit 1
fi

if [ -z "$2" ] ; then
	echo "Didn't gie in put Batch"
	exit 1
fi
if (( $2 <= 0 || $2 >= 11 )) ; then
	echo "Batch $2 isn't in [1,10]"
	exit 1
fi

job_name="$1"_batch_"$2"
vcf_file=vcfs/"$job_name".vcf.gz
date=$(date '+%Y%m%d_%H%M%S')
if [ -f "$vcf_file" ] ; then
	mv "$vcf_file" failed/"$job_name"_"$date".vcf.gz
fi
for file in output/"$job_name".* ; do 
	mv "$file" "failed/${job_name}_${date}.$(echo "$file" | cut -f2 -d/ | cut -f2- -d.)"
done

echo relaunching "$1" "$2"

var_str="STR=$1,BATCH=$2"
if [ -n "$3" ] ; then
	var_str="$var_str",CUT_SAMPLES="$3"
fi
if [ -n "$MAX_SAMPLES" ] ; then
	var_str="$var_str",MAX_SAMPLES="$MAX_SAMPLES"
fi
qsub -v "$var_str" call_strs_big.pbs > /dev/null
