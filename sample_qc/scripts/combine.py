import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-r', action='append')
parser.add_argument('-k', action='append')
parser.add_argument('-n', required = True)

args = parser.parse_args()

sample_dict = {}

first_keep = True
for file_name in args.k:
	with open(os.environ['UKB'] + '/sample_qc/snp_informed/hap/keep/' + file_name + '.sample') as keep_file:
		header = True
		current_keep = {}
		for line in keep_file:
			id = line.split()[0]
			if header :
				if id != "ID_1":
					print("Expected first line to be a header line in file " + file_name \
						+ " instead see " + id)
					exit(-1)
				header = False
				continue
			current_keep[id] = line
		if first_keep:
			sample_dict = current_keep
			first_keep = False
		else:
			for id in keep_dict:
				if id not in current_keep:
					del sample_dict[id]

for file_name in args.r:
	with open(os.environ['UKB'] + '/sample_qc/snp_informed/hap/remove/' + file_name + '.sample') as remove_file:
		header = True
		for line in remove_file:
			id = line.split()[0]
			if header:
				if id != "ID_1":
					print("Expected first line to be a header line in file " + file_name \
						+ " instead see " + id)
					exit(-1)
				header = False
				continue
			if id in sample_dict:
				del sample_dict[id]

with open(os.environ['UKB'] + '/sample_qc/snp_informed/hap/combined/' + args.n + '.sample', 'w') as outfile:
	outfile.write('ID_1 ID_2 missing sex\n')
	for id in sample_dict:
		outfile.write(sample_dict[id])
				
			
