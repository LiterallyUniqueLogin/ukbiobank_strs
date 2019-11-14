import os 
import argparse
import datetime
import logging
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument('sample_location', 
		help='name of the file (w/o prefix) in ./combined to verify')
args = parser.parse_args()

all_samples = {}
with open(os.environ['UKB'] + '/sample_qc/snp_informed/hap/combined/' + args.sample_location + '.sample') as sample_file:
	header = True
	for line in sample_file:
		id = line.split()[0]
		if header :
			if id != "ID_1":
				print("Expected first line to be a header line in file " + sample_file \
					+ " instead see " + id)
				exit(-1)
			header = False 
			continue
		all_samples[id] = line


with open(os.environ['UKB'] + '/non_genetic_data/ukbgene/ukb46122_rel_s488282.dat') as kinship_file:
	first = True
	for line in original_kinship_file:
		if first:
			first=False
			continue
		
		sample1, sample2 = line.split()[0:2]
		if sample1 in all_samples and sample2 in all_samples:
			print("Found related pair in the file", sample1, sample2)
			exit(-1)

print("No relatedness found")

