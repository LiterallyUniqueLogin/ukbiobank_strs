#!/bin/bash

if [ -z "$1" ] ; then
        echo "didn't give first argument - should be the run name"
        exit -1
fi

if [ -z "$2" ] ; then
        echo "didn't give second argument - should be chrom number"
        exit -1
fi

if [ -z "$3" ] ; then
        echo "didn't give third argument - should be number of lines to skip from the top of \
the callset vcf till you get to the first variant line"
        exit -1
fi

if [ -z "$TMPDIR" ] ; then
        echo "Didn't set TMPDIR"
        exit -1
fi

sed -e "s/%RUN_NAME%/$1/g" \
	$UKB/str_imputed/check_beagle_output.pbs \
	> $TMPDIR/check_beagle_output_$1.pbs

qsub -v "INPUT1=$2,INPUT2=$3" $TMPDIR/check_beagle_output_$1.pbs

