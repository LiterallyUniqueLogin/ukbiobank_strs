import os
import argparse
import subprocess as sp
import tempfile
import matplotlib.pyplot as plot;
import numpy as np
import time

#TODO
#use arg parse to figure out if we're testing imputation with full ref panel
#or with just eur samples
eur_only = True

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
	file_dir = ukb + '/pre_imputation_qc/ref_panel_size/output/variant_overlaps/chr' + chrom
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

#globals for the function below
#indexed by chromosome
cached_hipstr_lines = {}
nloci = {}

#Returns a list of calls
#Each call is a pair (dist, idx)
#Dist is the distribution of how likely 
#beagle thought each genotype is
#idx is the index of which genotype is 
#correct according to hipstr
#or -1 if the index is none of those
def gather_data(chrom, positions, sample):
	with tempfile.NamedTemporaryFile(mode = 'w+') as pos_file:
		positions_file_contents = chrom + "\t" + "\n1\t".join(positions) + "\n"
		pos_file.write(positions_file_contents)
		#with open("/home/jmargoli/temp.positions", 'w') as temp_pos:
		#	temp_pos.write(positions)
		positions = set(positions)

		if eur_only:
			panel_name = 'eur_panel'
		else:
			panel_name = 'full_panel'

		bcftools_command = "source ~/.bashrc && conda activate bcftools && "

		imputed_file_loc = ukb + '/pre_imputation_qc/ref_panel_size/{}_imputed/chr{}/{}.vcf.gz'.format(panel_name, chrom, sample)
		imputed_command = 'bcftools query -R ' + pos_file.name + ' -f "%POS %ID %REF %ALT [%GT %AP1 %AP2 ]\n" ' + imputed_file_loc
		imputed_lines = iter(sp.run(bcftools_command + imputed_command, shell = True, stdout = sp.PIPE, stderr = sp.PIPE, universal_newlines = True).stdout.split('\n'))

		#only one call to bcftools query is need for the entire chromosome,
		#so cached the results
		global cached_hipstr_lines, nloci
		if chrom not in cached_hipstr_lines:
			hipstr_file_loc = os.environ['DATASETS'] + '/1000Genomes/hipstr_calls_sample_subset/eur/chr{}/hipstr.chr{}.eur.vcf.gz'.format(chrom, chrom)
			hipstr_command = 'bcftools query -R ' + pos_file.name + ' -f "%POS %ID %REF %ALT [%SAMPLE %GT ]\n" ' + hipstr_file_loc
			cached_hipstr_lines[chrom] = sp.run(bcftools_command + hipstr_command, shell = True, stdout = sp.PIPE, stderr = sp.PIPE, universal_newlines = True).stdout.split('\n')
			nloci[chrom] = len(cached_hipstr_lines[chrom])

		#print(cached_hipstr_lines[chrom][:10])
		#exit()
		my_sample_idx = cached_hipstr_lines[chrom][0].split().index(sample)
		def hipstr_lines_iter():
			for line in cached_hipstr_lines[chrom]:
				line = line.split()
				if len(line) < 1:
					return None
				#filter to only our sample
				yield [line[idx] for idx in [0,1,2,3,my_sample_idx+1]]
		
		hipstr_lines = hipstr_lines_iter()

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
			correct_alleles = set() 
			for hip_allele_idx in map(int, hipstr_line[4].split('|')):
				if hipstr_alleles[hip_allele_idx] in imputed_alleles:
					correct_alleles.add(imputed_alleles.index(hipstr_alleles[hip_allele_idx]))
				else:
					correct_alleles.add(-1)
			
			imputedAP1 = list(map(float, imputed_line[5].split(',')))
			imputedAP1.insert(0, 1-np.sum(imputedAP1))
			imputedAP2 = list(map(float, imputed_line[6].split(',')))
			imputedAP2.insert(0, 1-np.sum(imputedAP1))

			#POS, imputed phased GT, imputed AP1, imputed AP2, correct unphased GT
			#Correct unphased GT is a set which contains two elements if the sample is heterozygous
			#and one if the sample is homozygous
			yield imputed_line[0], list(map(int, imputed_line[4].split('|'))), imputedAP1, imputedAP2, correct_alleles


#main method
sample_list = get_sample_list()

#these are totaled across all samples and chromosomes
ncalls = 0
ncorrect = 0
for chrom in range(1,23):
	chrom = str(chrom)
	positions = find_overlapping_STRs(chrom)

	loci_acc = {}
	for position in positions:
		loci_acc[position] = (0,0)

	for sample_num, sample in enumerate(sample_list):
		print("Processing chrom {} sample {} ({})".format(chrom, sample_num + 1, sample))
		for pos, imp_phased_gt, imp_ap1, imp_ap2, correct_unphased_gt in gather_data(chrom, positions, sample):
			ncalls += 1

			unphased_imp_call = {np.argmax(imp_ap1), np.argmax(imp_ap2)}
			if unphased_imp_call == correct_unphased_gt:
				ncorrect += 1
		break

	#for each ppt, gather those loci
	#analyze the results
	break

print("Number of overlapping loci", np.sum(list(nloci.values())))
print("Number of total calls (all samples, all loci)", ncalls)
print("Fraction of all (loci, sample) pairs hipstr declined to call {:.2%}".format(1 - ncalls/(49*np.sum(list(nloci.values())))))
print("Fraction of hipstr calls that were correctly imputed {:.2%}".format(ncorrect/ncalls))

#make the plots
#plots are histograms of predicted accuracy vs actual accuracy
#with error bars denoting what the expected value is
#one plot with just the mostly likely predicted values (which are not phased)
#one plot comparing calls (which are phased)
#one plot with all predicted value
#the above plots including and excluding calls that are not in the reference
#one plot with a distribution of the accuracy of different loci

#also need to do some comparison of how many STRs we're leaving out
