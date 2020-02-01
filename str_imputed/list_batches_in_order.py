import glob
import os
import argparse
import functools

#files must already exist to be listed

ukb = os.environ['UKB']

def list_of_files(run_name, chr):
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
	return sorted(files, key=functools.cmp_to_key(num_string_comparator))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("run_name", help="will look for the vcfs to merge in {ukb}/str_imputed/{run_name}/vcf_batches/chr{chr}_samples_*.vcf.gz")
	parser.add_argument("chr", help="the chromosome number")

	args=parser.parse_args()
	
	files = list_of_files(args.run_name, args.chr)

	for file in files:
		print(file)
