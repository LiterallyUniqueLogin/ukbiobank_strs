import gzip
import os
import sys
import glob
import argparse
import sys

ukb = os.environ['UKB']

def check_samples(run_name, chrom, printing=True):
	samples = ""
	with open(glob.glob(f"{ukb}/microarray/*.sample")[0], 'rt') as samples_file:
		samples_iter = iter(samples_file)
		next(samples_iter)
		first = True
		for line in samples_iter:
			line_head = line.split()[0]
			if line_head == "0":
				continue
			if first:
				first = False
			else:
				samples += "\t"
			samples += f"{line_head}_{line_head}"

	with gzip.open(f"{ukb}/str_imputed/runs/{run_name}/vcfs/chr{chrom}.vcf.gz", 'rt') as chrom_vcf:
		chrom_iter = iter(chrom_vcf)
		line = '##'
		while line[:2] == '##':
			line = next(chrom_iter)
		chrom_samples = line[46:].strip()
		if samples != chrom_samples :
			if not printing:
				return False
			print(f"Samples in chromosome {chrom} not the same as original callset samples!", file=sys.stderr)
			print("Original samples: ", samples, file=sys.stderr)
			print("Output samples:", chrom_samples, file=sys.stderr)
			exit(-1)
	
	if not printing:
		return True
	print(f"Success! Samples are correct for the final output chromosome {chrom} file")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("run_name", help="The name of this run of imputation")
	parser.add_argument("chromosome_number", help="The chromosome number to check")
	args = parser.parse_args()

	check_samples(args.run_name, args.chromosome_number)	

