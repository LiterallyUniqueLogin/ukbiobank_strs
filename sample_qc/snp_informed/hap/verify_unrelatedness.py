import os 
import argparse
import datetime
import logging
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument('sample_location', help='name of the file (w/o prefix or _unrelated) in ./combined to verify')
args = parser.parse_args()

def load_samples(file_loc):
	samples = {}
	with open(file_loc) as sample_file:
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
			samples[id] = line
	return samples

unfiltered_samples = load_samples(os.environ['UKB'] + '/sample_qc/snp_informed/hap/combined/' + args.sample_location +           '.sample')
filtered_samples   = load_samples(os.environ['UKB'] + '/sample_qc/snp_informed/hap/combined/' + args.sample_location + '_unrelated.sample')

neighbors = {}

for sample in unfiltered_samples:
	neighbors[sample] = set()

#Check if there is any relatedness
with open(os.environ['UKB'] + '/non_genetic_data/ukbgene/ukb46122_rel_s488282.dat') as kinship_file:
	first = True
	for line in kinship_file:
		if first:
			first=False
			continue
		
		sample1, sample2 = line.split()[0:2]
		if sample1 in filtered_samples and sample2 in filtered_samples:
			print("Found related pair in the file", sample1, sample2)
			exit(-1)
		if sample1 in unfiltered_samples and sample2 in unfiltered_samples:
			neighbors[sample1].add(sample2)
			neighbors[sample2].add(sample1)

#Check for maximality
for sample in unfiltered_samples:
	if sample in filtered_samples:
		continue
	error = True
	for neighbor in neighbors[sample]:
		if neighbor in filtered_samples:
			error = False
			break
	if error:
		print("Found a sample that was filtered but is not related to anyone left", sample)
		exit(-1)

#check for subset
for sample in filtered_samples:
	if sample not in unfiltered_samples:
		print("Found a sample in the unrelated set that's not in the original set", sample)
		exit(-1)

print("Success! No relatedness found. Remaining unrelated sample set is a maximal subset of the original.")


