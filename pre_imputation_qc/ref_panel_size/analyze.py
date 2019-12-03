import os
import argparse
import subprocess as sp
import tempfile

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

#Returns a list of calls
#Each call is a pair (dist, idx)
#Dist is the distribution of how likely 
#beagle thought each genotype is
#idx is the index of which genotype is 
#correct according to hipstr
#or -1 if the index is none of those
def gather_data(chrom, positions, sample):
	with tempfile.NamedTemporaryFile(mode = 'w+') as pos_file:
		positions = chrom + "\t" + "\n1\t".join(positions) + "\n"
		pos_file.write(positions)
		
		if eur_only:
			panel_name = 'eur_panel'
		else:
			panel_name = 'full_panel'

		bcftools_command = "source ~/.bashrc && conda activate bcftools && "

		imputed_file_loc = ukb + '/pre_imputation_qc/ref_panel_size/{}_imputed/chr{}/{}.vcf.gz'.format(panel_name, chrom, sample)
		imputed_command = 'bcftools query -R ' + pos_file.name + ' -f "%POS %REF %ALT [%GT %AP1 %AP2]\n" ' + imputed_file_loc
		imputed_lines = sp.run(bcftools_command + imputed_command, shell = True, stdout = sp.PIPE, stderr = sp.PIPE, universal_newlines = True).stdout
		print(imputed_lines[0:10])
		exit()

		hipstr_file_loc = os.environ['DATASETS'] + '/1000G/hipstr_calls_sample_subset/eur/chr{}/hipstr.chr{}.eur.vcf.gz'.format(chrom, chrom)
		hipstr_command = 'bcftools query -R ' + pos_file.name + ' -f "%POS %REF %ALT [%GT]\n" -s ' + sample + ' ' + hipstr_file_loc
		hipstr_lines = sp.run(bcftools_command + hipstr_command, shell = True, stdout = sp.PIPE, stderr = sp.PIPE, universal_newlines = True).stdout

		calls = []

		while True:
			imputed_line = next(imputed_lines)
			hipstr_line = next(hipstr_lines)

			if hipstr_line is None and imputed_line is not None:
				print("more lines in imputed file than hipstr file!")
				exit(-1)
			if hipstr_line is not None and imputed_line is None:
				print("More lines in hipstr file than imputed file!")
				exit(-1)
			if hipstr_line is None and imputed_line is None:
				break

			hipstr_line = hipstr_line.split()
			imputed_line = imputed_line.split()
	
			if hipstr_line[0] != imputed_line[0]:
				print("Different positions for hipstr and imputed lines! {} {}".format(hipstr_line[0], imputed_line[0]))
				exit(-1)

			yield imputed_line, hipstr_line


sample_list = get_sample_list()

for chrom in range(1,23):
	chrom = str(chrom)
	positions = find_overlapping_STRs(chrom)

	loci_acc = {}
	for position in positions:
		loci_acc[position] = (0,0)

	for sample in sample_list:
		for imputed_call, hipstr_call in gather_data(chrom, positions, sample):
			
	exit()

	#for each ppt, gather those loci
	#analyze the results
	break

#make the plots
#plots are histograms of predicted accuracy vs actual accuracy
#with error bars denoting what the expected value is
#one plot with just the mostly likely predicted values (which are not phased)
#one plot comparing calls (which are phased)
#one plot with all predicted value
#the above plots including and excluding calls that are not in the reference
#one plot with a distribution of the accuracy of different loci

#also need to do some comparison of how many STRs we're leaving out
