import os
import subprocess as sp
import argparse
import gzip 
import sys

parser = argparse.ArgumentParser()
parser.add_argument("vcf1", help="Location of the first vcf.gz file")
parser.add_argument("vcf2", help="Location of the second vcf.gz file")
parser.add_argument("--skip_lines1", help="Number of lines to skip at the beginning of the first vcf file")
parser.add_argument("--skip_lines2", help="Number of lines to skip at the beginning of the second vcf file")

args=parser.parse_args()

with gzip.open(args.vcf1, "rt") as vcf1, gzip.open(args.vcf2, "rt") as vcf2:

	if args.skip_lines1 is not None:
		for _ in range(int(args.skip_lines1)):
			next(vcf1)

	if args.skip_lines2 is not None:
		for _ in range(int(args.skip_lines2)):
			next(vcf2)

	count = 0
	while True:
		line1 = next(vcf1, None)
		line2 = next(vcf2, None)

		if line1 is None and line2 is None:
			print("Done. Files are identical")
			exit(0)
			
		if line1 is None and line2 is not None:
			print(f"Difference: vcf1 had fewer lines than vcf2, vcf1 terminated at line {count}", file=sys.stderr)
			exit(-1)
		
		if line2 is None and line1 is not None:
			print(f"Difference: vcf2 had fewer lines than vcf1, vcf2 terminated at line {count}", file=sys.stderr)
			exit(-1)
		
		count += 1

		line1_del = False
		line2_del = False

		line1 = line1.split()
		if len(line1) > 7:
			del line1[7]
			line1_del = True

		line2 = line2.split()
		if len(line2) > 7:
			del line2[7]
			line2_del = True

		if line1_del != line2_del:
			print(f"Difference: vcfs differ at line {count}, lines have spacing. vcf1 line: '{line1}', vcf2 line: '{line2}'", file=sys.stderr)
			exit(-1)

		if line1 != line2:
			print(f"Difference: vcfs differ at line {count}. vcf1 line: '{line1}', vcf2 line: '{line2}'", file=sys.stderr)
			exit(-1)
	
		if count % 1000 == 0:
			print(f"No difference in the first {count} lines", end="\r")
			sys.stdout.flush()

