
if [ -z "$BATCH" ] ; then echo "BATCH variable is unset" ; exit 1 ; fi
if [ -z "$STR" ] ; then echo "STR variable is unset" ; exit 1 ; fi
if (( ! ("$BATCH" >= 1 && "$BATCH" <= 10) )) ; then
	echo "BATCH $BATCH is not between 1 and 10" ; exit 1
fi

CLEANED_STR=$(echo "$STR" | sed -e 's%/%_%g')
STR_BED="$UKB"/exome/fe_cram_str_calls/str_beds/"$CLEANED_STR".bed
if [ ! -f "$STR_BED"  ] ; then
        echo "STR $STR (cleaned $CLEANED_STR) does not corerspond to a bed file in $UKB/exome/fe_cram_str_calls/str_beds/." 
	exit 1
fi

vcf_name=vcfs/"$STR"_batch_"$BATCH".vcf.gz

# file doesn't exist
if [ ! -f "$vcf_name" ] ; then
	echo "$STR"_batch_"$BATCH" file doesn\'t exist
	exit 1
fi

# file won't unzip properly
set -o pipefail
last_line=$(zcat "$vcf_name" 2> /dev/null | tail -n 1)
if (( $? != 0 )) ; then
	echo "$STR"_batch_"$BATCH" unzip fail
	exit 1
fi

# file truncated
if [[ "#" == $(echo "$last_line" | cut -c1-1) ]] ; then
	echo "$STR"_batch_"$BATCH" file truncated
	exit 1
fi

samples_per_batch=4000
# last line truncated
if (( $(echo "$last_line" | wc -w) != $(zcat "$vcf_name" | tail -n 2 | head -n 1 | wc -w) )) ; then
	echo "$STR"_batch_"$BATCH" last line truncated
	exit 1
fi

