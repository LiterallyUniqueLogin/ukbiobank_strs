#!/bin/bash

if [ -z "$1" ] ; then
        echo "didn't give first argument - should be the run name" 1>&2
        exit -1
fi

if [ -z "$2" ] ; then
        echo "didn't give second argument - should be chrom number" 1>&2
        exit -1
fi

if [ -z "$TMPDIR" ] ; then
        echo "Didn't set TMPDIR" 1>&2
        exit -1
fi

sed -e "s/%RUN_NAME%/$1/g" -e "s/%CHROM%/$2/g" \
	$UKB/str_imputed/scripts/concat_chr_regions.pbs \
	> $TMPDIR/concat_chr_regions_$1_$2.pbs

qsub $TMPDIR/concat_chr_regions_${1}_${2}.pbs

