#!/bin/bash

now=$(date +%Y_%M_%d_%H_%M_%S)
outdir=$UKB/sample_qc/primus_test/output_$now
mkdir outdir
/projects/ps-gymreklab/resources/source/PRIMUS_v1.9.0/bin/run_PRIMUS.pl \
	-i FILE=$UKB/non_genetic_data/ukbgene/ukb46122_rel_s488282.dat \
	FID1=1 IID1=1 FID2=2 IID2=2 PI_HAT=5 \
	--no_PR -t 0.08838835 -o outdir \
	> $outdir/stdout.txt 2> $outdir/stderr.txt
