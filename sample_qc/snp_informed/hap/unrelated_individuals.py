import os 
import argparse
import datetime
import logging
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument('sample_location', 
		help='name of the file (w/o prefix) in ./combined to take the unrelated sample subset of')
args = parser.parse_args()

now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
output_dir = os.environ['UKB'] + '/sample_qc/snp_informed/hap/primus_output/' + now
os.makedirs(output_dir, exist_ok = True)
logging.basicConfig(filename = output_dir + "/script_log.log", level=logging.DEBUG)

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
logging.info("Done reading original sample file. Num samples: {}".format(len(all_samples)))

#all sample numbers that appear in the related samples file
related_samples = set()

kinship_subset_filename = output_dir + "/kinship_subset.dat"
with open(os.environ['UKB'] + '/non_genetic_data/ukbgene/ukb46122_rel_s488282.dat') as original_kinship_file:
	with open(kinship_subset_filename, 'w') as kinship_subset_file:
		first = True
		for line in original_kinship_file:
			if first:
				first=False
				kinship_subset_file.write(line)
				continue
			
			sample1, sample2 = line.split()[0:2]
			if not (sample1 in all_samples and sample2 in all_samples):
				continue
			related_samples.add(sample1)
			related_samples.add(sample2)
			kinship_subset_file.write(line)

logging.info("Done creating subset kinship file. Num related samples: {}".format(len(related_samples)))

commandString = \
os.environ['SOURCE'] + "/PRIMUS_v1.9.0/bin/run_PRIMUS.pl " + \
"-i FILE=$UKB/non_genetic_data/ukbgene/ukb46122_rel_s488282.dat " + \
"FID1=1 IID1=1 FID2=2 IID2=2 PI_HAT=5 " + \
"--no_PR -t 0.08838835 -o " + output_dir+"_PRIMUS"
#the -t threshold as described in the UKB nature paper

output = sp.run(commandString, shell = True, stdout = sp.PIPE, stderr = sp.PIPE)
logging.info("Done running PRIMUS")
logging.info("PRIMUS output: " + output.stdout.decode())
logging.info("PRIMUS error (if any): " + output.stderr.decode())

with open(os.environ['UKB'] + '/sample_qc/snp_informed/hap/combined/' + args.sample_location + '_unrelated.sample', 'w') as output:
	output.write("ID_1 ID_2 missing sex\n")
	#write out all the samples that aren't related to anyone
	for sample in all_samples:
		if sample not in related_samples:
			output.write(all_samples[sample])

	#write out all the samples that were selected among the related ones	
	with open(output_dir + "_PRIMUS/ukb46122_rel_s488282.dat_maximum_independent_set") as unrelated_file:
		first = True
		for line in unrelated_file:
			if first:
				first = False
				continue
			output.write(all_samples[line.strip()])

