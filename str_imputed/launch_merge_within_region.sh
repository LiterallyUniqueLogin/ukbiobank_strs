#!/bin/bash

if [ -z "$1" ] ; then
        echo "didn't give first argument - should be the run name"
        exit -1
fi

if [ -z "$2" ] ; then
        echo "didn't give second argument - should be chrom number"
        exit -1
fi

if [ -z "$TMPDIR" ] ; then
        echo "Didn't set TMPDIR"
        exit -1
fi

sed -e "s/%RUN_NAME%/$1/g" \
	$UKB/str_imputed/merge_within_region.pbs \
	> $TMPDIR/merge_within_region_$1.pbs

for iter in $(seq 0 49) ; do
	start_pos=$(( iter*5000000 + 1 ))
	qsub -v "INPUT1=${start_pos},INPUT2=$2" $TMPDIR/merge_within_region_$1.pbs
done

