import glob
import os
import argparse
import functools

parser = argparse.ArgumentParser()
parser.add_argument("run_name", help="will look for the vcfs to merge in {ukb}/str_imputed/{run_name}/vcf_batches/chr{chr}_samples_*.vcf.gz")
parser.add_argument("chr", help="the chromosome number")

args=parser.parse_args()
chr = args.chr
run_name = args.run_name

ukb = os.environ['UKB']

files = glob.glob(f"{ukb}/str_imputed/{run_name}/batches/chr{chr}_samples_*.vcf.gz")
def num_string_comparator(a, b):
	if len(a) != len(b):
		return len(a) - len(b)
	elif a<b:
		return -1
	elif a == b:
		return 0
	else:
		return 1
files = sorted(files, key=functools.cmp_to_key(num_string_comparator))

for file in files:
	print(file)

