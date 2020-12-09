#!/usr/bin/env bash

CLEANED_STR=$(echo "$STR" | sed -e 's%/%_%g')

if [ -z "$BATCH" ] ; then echo "BATCH variable is unset" ; exit 1 ; fi
if [ -z "$STR" ] ; then echo "STR variable is unset" ; exit 1 ; fi
if (( ! ("$BATCH" >= 1 && "$BATCH" <= 10) )) ; then
	echo "BATCH $BATCH is not between 1 and 10" ; exit 1 
fi

STR_BED="$UKB"/exome/fe_cram_str_calls/str_beds/"$CLEANED_STR".bed
if [ ! -f "$STR_BED"  ] ; then
	echo "STR $STR (cleaned $CLEANED_STR) does not corerspond to a bed file in $UKB/exome/fe_cram_str_calls/str_beds/." ; exit 1 
fi

ori_max_reads="$MAX_READS"
if [ -z "$MAX_READS" ] ; then
	MAX_READS=1000000
elif [[ ! $MAX_READS =~ [0-9]+ ]] ; then
	echo "MAX_READS $MAX_READS is not a positive integer"
	exit 1
fi

echo BATCH "$BATCH" STR "$STR" CUT_SAMPLES "$CUT_SAMPLES" MAX_READS "$ori_max_reads"
echo BATCH "$BATCH" STR "$STR" CUT_SAMPLES "$CUT_SAMPLES" MAX_READS "$ori_max_reads" 1>&2 

job_name="$STR"_batch_"$BATCH"

cd "$UKB"/exome/fe_crams/ || { echo "Can't get to $UKB/exome/fe_crams/ ; exiting " ; exit 1 ; }

# keep track of memory usage
{ while true; do
	mem_monitor
	sleep 5
	echo
done } > "$UKB"/exome/fe_cram_str_calls/output/"$job_name".mem_monitor &

samples_per_batch=4000
CRAM_FILES=$({ for file in *cram ; do if [[ ! "$file" =~ _ ]] ; then echo "$file" ; fi ; done ; } | head -n $((samples_per_batch*BATCH)) | tail -n "$samples_per_batch" | tr '\n' , | sed -e 's/,$//')
# cut samples should be a dash delimited list of sample names
for sample in $(echo "$CUT_SAMPLES" | sed -e 's/-/ /g'); do
	CRAM_FILES="$(echo "$CRAM_FILES" | sed -e 's/,'"$sample"'.cram,/,/' -e 's/^'"$sample"'.cram,//' -e 's/,'"$sample"'.cram$//')"
done
echo CRAM_FILES "$CRAM_FILES"

#uses hg38
#https://biobank.ctsu.ox.ac.uk/showcase/refer.cgi?id=1000
#need to set the library for each cram file because the cram file RG tags don't have LB fields
#so just use the CRAM filename for the library name.
# only 4000 files allowed open at once
# see /etc/security/limits.conf
{ echo BATCH "$BATCH" PBS_JOBID "$PBS_JOBID" CUT_SAMPLES "$CUT_SAMPLES" MAX_READS "$ori_max_reads"
  echo BATCH "$BATCH" PBS_JOBID "$PBS_JOBID" CUT_SAMPLES "$CUT_SAMPLES" MAX_READS "$ori_max_reads" 1>&2
  time /projects/ps-gymreklab/jmargoli/ukbiobank/utilities/hipstr/hipstr/HipSTR \
	--bams "$CRAM_FILES" \
	--bam-libs "$CRAM_FILES" \
	--bam-samps "$CRAM_FILES" \
	--fasta "$UKB"/exome/fe_cram_str_calls/hg38_no_chr.fa \
	--regions "$STR_BED" \
	--str-vcf "$UKB"/exome/fe_cram_str_calls/vcfs/"$job_name".vcf.gz \
	--max-reads "$MAX_READS" ; } \
	> "$UKB"/exome/fe_cram_str_calls/output/"$job_name".hipstr.stdout \
	2> "$UKB"/exome/fe_cram_str_calls/output/"$job_name".hipstr.stderr

# kill the mem_monitor process
kill -9 %1

