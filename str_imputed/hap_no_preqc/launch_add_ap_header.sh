#!/bin/bash

if [ -z "$1" ] ; then
	echo "Expecting the chromosome number as the first argument"
	exit -1
fi

for begin in $(seq 2001 1000 487409); do 
	qsub -v "INPUT1=$begin,INPUT2=$1" add_ap_header.pbs
done
