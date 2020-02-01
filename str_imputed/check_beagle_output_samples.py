import os
import sys
import argparse
import list_batches_in_order
import gzip

#doesn't allow IDs that look like "0"
#doesn't allow underscores in sample names

parser = argparse.ArgumentParser()
parser.add_argument("run_name", help="The name of this run of imputation")
parser.add_argument("sample_file", help="the .sample file with \
the list of samples in this imputation run, see $UKB/microarray/*.sample for an example")
parser.add_argument("chromosome_number", help="The number of the chromosome you wish to check")

args = parser.parse_args()

vcf_files = list_batches_in_order.list_of_files(args.run_name, args.chromosome_number)

def except_bad_values_iter(sample_file):
	my_iter = iter(sample_file)
	while True:
		val = next(my_iter, None)
		if val is None:
			break
	
		if val.strip() == "":
			continue 

		val = val.split()[0]
		if val == "0":
			continue

		yield val

with open(args.sample_file) as samples_file:
	next(samples_file) #skip the header line
	samples_file = except_bad_values_iter(samples_file)
	for file_num, vcf_file in enumerate(vcf_files):
		file_num += 1
		print(f"Working with file {vcf_file}                 ", end="\r")
		with gzip.open(vcf_file, 'rt') as vcf:
			for _ in range(10):
				next(vcf)
			vcf_line = next(vcf)
			vcf_samples = iter(vcf_line.split()[9:])
			while True:
				vcf_sample = next(vcf_samples, None)
				if vcf_sample is None:
					break #go to next file
				vcf_sample = vcf_sample.split("_")[0]

				true_sample = next(samples_file, None)
				if true_sample is None:
					print(f"Found more samples in the vcfs than in the samples file.\
 Extra sample {vcf_sample} in file {vcf_file} (file num {file_num})", file = sys.stderr)
					exit(-1)
					
				
				if vcf_sample == true_sample:
					continue
				else:
					print(f"Found differing samples in a beagle output vcf\
 than in the samples file. Sample file sample {true_sample}, vcf sample {vcf_sample} in file\
 {vcf_file} (file num {file_num})", file = sys.stderr)
					exit(-1)

	samples_line = next(samples_file, None)
	if samples_line is not None and samples_line.strip() != "":
		print(f"Found more samples in samples file than in all the vcfs. \
Extra sample {samples_line.split()[0]} in samples file", file = sys.stderr)
		exit(-1)
		
	print("Success! All samples are represented in the same order in the beagle\
output files as in the original sample file.                                    ")
