import argparse
import os
import sys
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument("run_name", help="Output files will be put in $UKB/str_imputed/runs/run_name/batches")
parser.add_argument("chromosome_number", help="the number of the chromosome to regional merge")
parser.add_argument("--regions", nargs="+", help="a list of regions to rerun (no other regions will be run), e.g. --nargs 5000001 210000001 ")

args = parser.parse_args()
run_name = args.run_name
chrom = args.chromosome_number

def error(msg):
	print(msg, file=sys.stderr)
	exit(-1)

if not "TMPDIR" in os.environ:
	error("Didn't set TMPDIR environment variable")

tmpdir = os.environ["TMPDIR"]
ukb = os.environ["UKB"]

if args.regions:
	for region in args.regions:
		if (int(region) - 1) % 5000000 != 0:
			error("Each region should be an index that is some multiple of 5 million plus 1")

sp.run(f'sed -e "s/%RUN_NAME%/{run_name}/g" \
	{ukb}/str_imputed/merge_within_region.pbs \
	> {tmpdir}/merge_within_region_{run_name}.pbs', shell=True)

if args.regions:
	regions = args.regions
else:
	regions = range(1,250000001,5000000)

for start_pos in regions:
	sp.run(f'qsub -v "INPUT1={start_pos},INPUT2={chrom}" {tmpdir}/merge_within_region_{run_name}.pbs', shell=True)

