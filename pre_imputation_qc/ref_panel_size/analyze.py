import os
import argparse
import subprocess as sp
import tempfile
import matplotlib.pyplot as plot
import numpy as np
import time
import argparse
import pickle
from collections import OrderedDict

ukb = os.environ['UKB']

def get_sample_list():
	with open(ukb + "/pre_imputation_qc/ref_panel_size/output/compare_samples/eur_hipstr_and_ref.sample") as samples_file:
		next(samples_file)
		samples = []
		for line in samples_file:
			samples.append(line.strip())
		return samples

#returns a list of positions to look at
#will not consider varaints which are located at different
#positions but are alternate spellings of one another
#e.g. 10 AAA AAAA vs 12 A AA
#Also will not consider varaints with different ids between
#the two vcfs, or different ref alleles
def find_overlapping_STRs(chrom):
	file_dir = ukb + '/pre_imputation_qc/ref_panel_size/output/variant_pos_overlaps/chr' + chrom
	with open(file_dir + '/0000.vcf') as ref_panel_file, \
		open(file_dir + '/0001.vcf') as hipstr_file:

		#skip metadata lines
		line = next(ref_panel_file)
		while line[0:2] == '##':
			line = next(ref_panel_file)

		line = next(hipstr_file)
		while line[0:2] == '##':
			line = next(hipstr_file)

		#figure out the columns we're interested in
		header = line.split()
		pos_idx = header.index('POS')
		ref_idx = header.index('REF')
		alt_idx = header.index('ALT')
		id_idx = header.index('ID')
	
		next(ref_panel_file)	
		next(hipstr_file)

		positions = []

		#loop over all overlapping positions
		while True:
			ref_panel_line = next(ref_panel_file, None)
			hipstr_line = next(hipstr_file, None)

			if hipstr_line is None and ref_panel_line is not None:
				print("more lines in ref_panel file than hipstr file!")
				exit(-1)
			if hipstr_line is not None and ref_panel_line is None:
				print("More lines in hipstr file than ref_panel file!")
				exit(-1)
			if hipstr_line is None and ref_panel_line is None:
				break

			hipstr_line = hipstr_line.split()
			ref_panel_line = ref_panel_line.split()

			if hipstr_line[pos_idx] != ref_panel_line[pos_idx]:
				print("Non-matching positions in ref_panel and hipstr files!")
				exit(-1)
			pos = hipstr_line[pos_idx]

			if hipstr_line[id_idx] != ref_panel_line[id_idx]:
				continue

			if hipstr_line[ref_idx] != ref_panel_line[ref_idx]:
				continue

			positions.append(pos)

		return positions

#globals for the two functions below
#indexed by chromosome
nloci = {}
ngenotypes = {} #number of genotypes, summed across all loci
cached_hipstr_iters = {}
cached_hipstr_lines = {}
#indexed by chrom and sample
cached_imputed_iters = {}

def load_data(chrom, positions, samples):
	global cached_hipstr_iters, cached_hipstr_lines, cached_imputed_iters, nloci, eur_only

	bcftools_command = "source ~/.bashrc && conda activate bcftools && "

	if eur_only:
		panel_name = 'eur_panel'
	else:
		panel_name = 'full_panel'

	with tempfile.NamedTemporaryFile(mode = 'w+', delete = False) as pos_file:
		positions_file_contents = chrom + "\t" + ("\n" + chrom + "\t").join(positions) + "\n"
		pos_file.write(positions_file_contents)
		
		hipstr_file_loc = os.environ['DATASETS'] + '/1000Genomes/hipstr_calls_sample_subset/eur/chr{}/hipstr.chr{}.eur.vcf.gz'.format(chrom, chrom)
		hipstr_command = 'bcftools query -R ' + pos_file.name + ' -f "%POS %ID %REF %ALT [%SAMPLE %GT ]\n" ' + hipstr_file_loc
		cached_hipstr_iters[chrom] = sp.Popen(bcftools_command + hipstr_command, shell = True, stdout = sp.PIPE, stderr = sp.PIPE, universal_newlines = True).stdout

		cached_imputed_iters[chrom] = {}
		for sample in samples:
			imputed_file_loc = ukb + '/pre_imputation_qc/ref_panel_size/{}_imputed/chr{}/{}.vcf.gz'.format(panel_name, chrom, sample)
			imputed_command = 'bcftools query -R ' + pos_file.name + ' -f "%POS %ID %REF %ALT [%GT %AP1 %AP2 ]\n" ' + imputed_file_loc
			cached_imputed_iters[chrom][sample] = sp.Popen(bcftools_command + imputed_command, shell = True, stdout = sp.PIPE, stderr = sp.PIPE, universal_newlines = True).stdout
		
		#turn hipstr iter into strings 
		cached_hipstr_lines[chrom] = cached_hipstr_iters[chrom].readlines()
		nloci[chrom] = len(cached_hipstr_lines[chrom])

#Returns a list of calls
#Each call is a pair (dist, idx)
#Dist is the distribution of how likely 
#beagle thought each genotype is
#idx is the index of which genotype is 
#correct according to hipstr
#or -1 if the index is none of those
def read_data(chrom, positions, sample, record_ngenotypes):
	global cached_hipstr_lines, cached_imputed_iters, nloci, ngenotypes
	positions = set(positions)

	if record_ngenotypes:
		ngenotypes[chrom] = 0

	my_sample_idx = cached_hipstr_lines[chrom][0].split().index(sample)
	def hipstr_lines_iter():
		for line in cached_hipstr_lines[chrom]:
			line = line.split()
			if len(line) < 1:
				return None
			#filter to only our sample
			yield [line[idx] for idx in [0,1,2,3,my_sample_idx+1]]
	
	hipstr_lines = hipstr_lines_iter()
	imputed_lines = cached_imputed_iters[chrom][sample]

	last_pos = None
	call_idx = 0
	while True:
		imputed_line = next(imputed_lines, None)
		hipstr_line = next(hipstr_lines, None)
		call_idx += 1

		if call_idx == nloci[chrom] and next(imputed_lines, None) is None:
			#the last line of bcftools query will be a blank one
			break

		print("Processing call {}/{}".format(call_idx, nloci[chrom]), end='\r')

		if imputed_line is None:
			if hipstr_line is None:
				break
			else:
				print("more hipstr lines than expected!")
				exit(-1)

		imputed_line = imputed_line.split()

		#handle the fact that sometimes the imputed file has multiple
		#variants at the exact same position and we need to ignore
		#the non-STR one that's also being returned by the bcftools query
		while last_pos is not None and imputed_line[0] == last_pos:
			imputed_line = next(imputed_lines)
			imputed_line = imputed_line.split()

		#handle the fact that in the imputation panel
		#there are some STRs which overlap one another and have the 
		#same ID but are at different loci, some of which
		#we didn't ask for
		if imputed_line[0] not in positions:
			imputed_line = next(imputed_lines)
			imputed_line = imputed_line.split()

		if imputed_line is not None and hipstr_line is None:
			print("more imputed lines than expected!")
			exit(-1)

		last_pos = hipstr_line[0]

		#more handling the multiple variants same locus issue
		#Ignore variants that line up but have different IDs
		while hipstr_line[1] != imputed_line[1]:
			imputed_line = next(imputed_lines)
			imputed_line = imputed_line.split()

		if hipstr_line[0] != imputed_line[0]:
			print("Different positions for hipstr and imputed lines! hipstr pos: {} imputed pos: {} linenum: {}".format(hipstr_line[0], imputed_line[0], call_idx))
			exit(-1)

		if hipstr_line[4] == '.':
			#Hipstr couldn't genotype this sample at this locus
			continue

		imputed_alleles = imputed_line[3].split(',')
		imputed_alleles.insert(0, imputed_line[2])
		hipstr_alleles = hipstr_line[3].split(',')
		hipstr_alleles.insert(0, hipstr_line[2])

		if record_ngenotypes:
			nalleles = len(imputed_alleles)
			ngenotypes[chrom] += nalleles * (nalleles + 1)/2

		correct_alleles = set() 
		for hip_allele_idx in map(int, hipstr_line[4].split('|')):
			if hipstr_alleles[hip_allele_idx] in imputed_alleles:
				correct_alleles.add(imputed_alleles.index(hipstr_alleles[hip_allele_idx]))
			else:
				correct_alleles.add(-1)
		
		imputedAP1 = list(map(float, imputed_line[5].split(',')))
		imputedAP1.insert(0, max(1-np.sum(imputedAP1), 0))
		imputedAP2 = list(map(float, imputed_line[6].split(',')))
		imputedAP2.insert(0, max(1-np.sum(imputedAP2), 0))

		#POS, imputed phased GT, imputed AP1, imputed AP2, correct unphased GT
		#Correct unphased GT is a set which contains two elements if the sample is heterozygous
		#and one if the sample is homozygous
		yield imputed_line[0], list(map(int, imputed_line[4].split('|'))), imputedAP1, imputedAP2, correct_alleles


class CallAggregator:

	def __init__(self, total_calls, bin_size = 0.01): #default bin size of 1%
		self.calls = np.full(total_calls, np.nan)
		self.call_successes = np.zeros(total_calls, dtype=np.bool) #this is all falses
		self.ncalls = 0
		self.call_limit = total_calls
		self.bin_size = bin_size

		#bins are [bin_size*n, bin_size*(n+1)) except for the first bin which is (0, bin_size) and the last bin which is [bin_size*n, 100)
		self.call_bins = OrderedDict()    #total number of calls which claimed this accuracy
		self.success_bins = OrderedDict() #total number of calls which claimed this accuracy and were correct
		bins = list(np.arange(0, 1, self.bin_size))
		bins.insert(0, 'No') #imputed probability is 0
		bins.append('Yes') #imputed probability is 1
		for bin in bins:
			self.call_bins[bin] = 0
			self.success_bins[bin] = 0

	#predicted fraction of being correct, and whether it was correct or not	
	def add(self, fraction, correct):
		if self.ncalls >= self.call_limit:
			self.calls = np.resize(self.calls, (self.call_limit*2))
			self.call_successes = np.resize(self.call_successes, (self.call_limit*2))
			self.call_limit *= 2

		#handle recording
		self.calls[self.ncalls] = fraction
		self.call_successes[self.ncalls] = correct

		#handle bins	
		if fraction == 0:
			idx = 'No'
		elif fraction == 1:
			idx = 'Yes'
		else:	
			idx = (fraction // self.bin_size)*self.bin_size
			
		self.call_bins[idx] += 1
		if correct:
			self.success_bins[idx] += 1
		
		self.ncalls += 1

pickle_file_loc = 'analyze_data_{}.pickle'
def load_pickle(eur_only):
	if eur_only:
		file_loc = pickle_file_loc.format('eur')
	else:
		file_loc = pickle_file_loc.format('full')
	with open(file_loc, 'rb') as f:
		return pickle.load(f)

#main method
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--eur_only', action="store_true")
	args = parser.parse_args()
	eur_only = bool(args.eur_only)
	if eur_only:
		pickle_suffix = 'eur'
	else:
		pickle_suffix = 'full'

	sample_list = get_sample_list()

	#accrue information per sample
	sample_acc = {} #unphased
	for sample in sample_list:
		sample_acc[sample] = [0,0] #[number calls, number successes]

	#accrue info per locus
	loci_acc = {} #unphased

	#these are totaled across all samples and chromosomes
	ncalls = 0
	ncorrect_unphased = 0
	ncorrect_phased = 0
	ncalls_not_in_ref = 0
	ncalls_predicted = 9716850 #from previous runs

	#record calibration stats about calls beagle makes
	unphased_call_calibration = CallAggregator(ncalls_predicted) #the most likely unphased genotype

	  #the most likely phased genotype
	  #note that because Hipstr doesn't output phased genotypes,
	  #we call the phased genotype call a|b correct if hipstr outputs either a/b or b/a
	phased_call_calibration = CallAggregator(ncalls_predicted) 

	unphased_genotype_calibration = CallAggregator(10)

	positions = {} #indexed by chrom
	for chrom in range(1,23):
		#load the positions and data for this chrom
		print("Loading positions for chrom {}".format(chrom), end='\r')
		chrom = str(chrom)
		positions[chrom] = find_overlapping_STRs(chrom)
		load_data(chrom, positions[chrom], sample_list)

		loci_acc[chrom] = {}
		for position in positions[chrom]:
			loci_acc[chrom][position] = [0,0] #[number calls, number successes]

		#read the data in for this chrom
		for sample_num, sample in enumerate(sample_list):
			print("Processing chrom {} sample {} ({})".format(chrom, sample_num + 1, sample))
			record_ngenotypes = sample_num == 0
			for pos, imp_phased_gt, imp_ap1, imp_ap2, correct_unphased_gt in read_data(chrom, positions[chrom], sample, record_ngenotypes):
				#record the call location	
				loci_acc[chrom][pos][0] += 1
				sample_acc[sample][0] += 1
				ncalls += 1

				#is the unphased call correct?	
				unphased_imp_perc = None
				unphased_imp_call = None
				#consider all genotypes a|b or b|a
				for allele1 in range(len(imp_ap1)):
					for allele2 in range(allele1):
						perc = imp_ap1[allele1]*imp_ap2[allele2] + \
							imp_ap1[allele2]*imp_ap2[allele1]
						genotype = {allele1, allele2}
						correct = genotype == correct_unphased_gt
						unphased_genotype_calibration.add(perc, correct)
						if unphased_imp_perc is None or unphased_imp_perc < perc:
							unphased_imp_perc = perc
							unphased_imp_call = genotype
							unphased_correct = correct

				#consider all genotypes a|a
				for allele in range(len(imp_ap1)):
					perc = imp_ap1[allele] * imp_ap2[allele]
					genotype = {allele}
					correct = genotype == correct_unphased_gt
					unphased_genotype_calibration.add(perc, correct)
					if unphased_imp_perc < perc:
						unphased_imp_perc = perc
						unphased_imp_call = genotype
						unphased_correct = correct

				#if not unphased_correct:
				#	print("Incorrect call chr: {} sample: {} locus: {} best gt {} imp_ap1 {} imp_ap2 {} perc {:.2f} true gt {} ".format(chr, sample, pos, unphased_imp_call, imp_ap1, imp_ap2, unphased_imp_perc, correct_unphased_gt))

				if unphased_correct:
					ncorrect_unphased += 1
					loci_acc[chrom][pos][1] += 1
					sample_acc[sample][1] += 1
				elif -1 in correct_unphased_gt:
					ncalls_not_in_ref += 1

				unphased_call_calibration.add(unphased_imp_perc, unphased_correct)
	
				#is the phased call correct?
				phased_imp_call = set(imp_phased_gt)
				phased_correct = phased_imp_call == correct_unphased_gt
				if phased_correct:
					ncorrect_phased += 1
				phased_imp_perc = imp_ap1[imp_phased_gt[0]]*imp_ap2[imp_phased_gt[1]]
				phased_call_calibration.add(phased_imp_perc, phased_correct)

		#save the data we've collected so far
		with open(pickle_file_loc.format(pickle_suffix), 'wb') as pickle_file:
			data = {"nloci" : nloci,
				"ncalls" : ncalls,
				"ngenotypes" : ngenotypes,
				"ncorrect_unphased" : ncorrect_unphased, 
				"ncorrect_phased" : ncorrect_phased, 
				"ncalls_not_in_ref" : ncalls_not_in_ref,
				"sample_acc" : sample_acc,
				"loci_acc" : loci_acc,
				"unphased_call_calibration" : unphased_call_calibration,
				"phased_call_calibration" : phased_call_calibration,
				"unphased_genotype_calibration": unphased_genotype_calibration,
				"max_chrom": chrom,
				"done": chrom == '22'}
			pickle.dump(data, pickle_file, pickle.HIGHEST_PROTOCOL)

	print("Number of overlapping loci", np.sum(list(nloci.values())))
	print("Number of total calls (all samples, all loci)", ncalls)
	print("Fraction of all (loci, sample) pairs hipstr declined to call {:.2%}".format(1 - ncalls/(49*np.sum(list(nloci.values())))))
	print("""Fraction of imputed unphased calls that are correct according to hipstr \
		(ignoring pairs hiptsr declined to call) {:.2%}""".format(ncorrect_unphased/ncalls))
	print("""Fraction of imputed phased calls that are correct according to hipstr {:.2%} \
		 (ignoring pairs hiptsr declined to call) \
		 (Note that because hipstr returns unphased calls, \
		  we consider a call a|b accurate if hipstr returns either a|b or b|a. \
		  So this overestimates the phased correctness) """.format(ncorrect_phased/ncalls))
	print("Fraction of hipstr calls with at least one allele with a genotype not in the reference panel {:.2%}".format(ncalls_not_in_ref/ncalls))

#make the plots
#plots are histograms of predicted accuracy vs actual accuracy
#with error bars denoting what the expected value is
#one plot with just the mostly likely predicted values (which are not phased)
#one plot comparing calls (which are phased)
#one plot with all predicted values
#the above plots including and excluding calls that are not in the reference
#one plot with a distribution of the accuracy of different loci
#and the number of hipstr no-calls at different loci
#also visualize this is in IGV
#and one with a distribution of accuracy across samples

#also need to do some comparison of how many STRs we're leaving out
#what are the backgrounds of these samples?

#Look for patterns in missingness of loci per sample?
