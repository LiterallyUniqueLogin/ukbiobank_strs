import argparse
import os
import sys
import subprocess as sp


def do_launch(run_name, chrom, regions, command_line = True):
	def error(msg):
		if command_line:
			print(msg, file=sys.stderr)
			exit(-1)
		else:
			return msg

	if not "TMPDIR" in os.environ:
		error("Didn't set TMPDIR environment variable")

	tmpdir = os.environ["TMPDIR"]
	ukb = os.environ["UKB"]

	if regions:
		for region in args.regions:
			if (int(region) - 1) % 5000000 != 0:
				error("Each region should be an index that is some multiple of 5 million plus 1")

	if not regions:
		regions = range(1,250000000,5000000)

	for start_pos in regions:
		sp.run(f'sed -e "s/%RUN_NAME%/{run_name}/g" -e "s/%CHROM%/{chrom}/g" -e "s/%POS%/{start_pos}" \
			{ukb}/str_imputed/scripts/merge_within_region.pbs \
			> {tmpdir}/merge_within_region_{run_name}_{chrom}_{pos}.pbs', shell=True)
		sp.run(f'qsub {tmpdir}/merge_within_region_{run_name}_{chrom}_{pos}.pbs', shell=True)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("run_name", help="Output files will be put in $UKB/str_imputed/runs/run_name/batches")
	parser.add_argument("chromosome_number", help="the number of the chromosome to regional merge")
	parser.add_argument("--regions", nargs="+", help="a list of regions to rerun (no other regions will be run), e.g. --nargs 5000001 210000001 ")

	args = parser.parse_args()
	run_name = args.run_name
	chrom = args.chromosome_number
	do_launch(run_name, chrom, args.regions)

