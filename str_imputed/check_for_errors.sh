#!/bin/sh

if [ -z "$1" ] ; then
        echo "Didn't give first argument - should be chrom number"
        exit -1
fi

grep -h -b100 -v INPUT $UKB/str_imputed/runs/$1/batches/output/*.e* | vim -

