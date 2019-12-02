import os
import argparse
import subprocess as sp

#TODO
#use arg parse to figure out if we're testing imputation with full ref panel
#or with just eur samples
eur_only = True

ukb = os.environ['UKB']

#returns a list of positions to look at
#will not consider varaints which are located at different
#positions but are alternate spellings of one another
#e.g. 10 AAA AAAA vs 12 A AA
def find_overlapping_STRs(chrom):
	file_dir = ukb + '/pre_imputation_qc/ref_panel_size/output/variant_overlap/chr' + chrom 
	with open(file_dir + '/0000.vcf') as ref_panel_file,
		open(file_dir + '/0001.vcf') as hipster_file:

		#skip metadata lines
		line = ref_panel_file.read_line()
		while line[0:2] == '##':
			line = ref_panel_file.read_line()

		line = hipstr_file.read_line()
		while line[0:2] == '##':
			line = hipstr_file.read_line()

		#figure out the columns we're interested in
		header_line = ref_panel_file.read_line()
		header = header_line.split()
		pos_idx = line.indexOf('POS')
		ref_idx = line.indexOf('REF')
		alt_idx = line.indexOf('ALT')
		
		hipstr_file.read_line()

		positions = []

		#loop over all overlapping positions
		while True:
			ref_panel_line = ref_panel_file.read_line()
			hipstr_line = hipstr_file.read_line()

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

			hipstr_alleles = set() 
			hipstr_alleles.add(hipstr_line[ref_idx])
			hipstr_alleles.add_all(hipstr_line[alt_idx].split(','))

			ref_panel_alleles = set() 
			ref_panel_alleles.add(ref_panel_line[ref_idx])
			ref_panel_alleles.add_all(ref_panel_line[alt_idx].split(','))

			if ref_panel_alleles.intersection(hipstr_alleles).is_empty():
				#no overlap, assume these are different variants
				#identified at the same position and ignore them
				continue

			hip_minus_ref = hipstr_alleles - ref_panel_alleles
			if not hip_minus_ref.empty():
				print("At pos {}, found alleles {} in hipstr that aren't represented in the ref panel".format(pos, hip_minus_ref))
				exit(-1)

			#Otherwise we know hipstr output is a subset of what is seens as possible in the ref panel)

			if not hipstr_line[ref_idx] == ref_panel_line[ref_idx]:
				print("Hipstr and ref_panel have overlapping alleles but different ref alleles at pos {}".format(pos))
				exit(-1)

			positions.add(pos)

		return positions

#Returns a list of calls
#Each call is a pair (dist, idx)
#Dist is the distribution of how likely 
#beagle thought each genotype is
#idx is the index of which genotype is 
#correct according to hipstr
#or -1 if the index is none of those
def gather_data(chrom, positions, sample):
	positions_arg = chrom + ':' + (',' + chrom  + ':').join(positions)
	
	if eur_only:
		panel_name = 'eur_panel'
	else:
		panel_name = 'full_panel'

	imputed_file_loc = ukb + '/pre_imputation_qc/ref_panel_size/()/chr{}/{}.vcf.gz'.format(full_panel, chrom, sample)
	imputed_command = 'bcftools query -r ' + positions_arg + ' -f "%REF %ALT [%GT %FORMAT/GP]\n" ' + imputed_file_loc  #TODO the format notation
	imputed_lines = sp.run(imputed_command, shell = True, stdout = sp.PIPE, stderr = sp.PIPE).stdout.decode().read_lines()

	hipstr_file_loc = os.environ['DATASETS'] + '/1000G/hipstr_only.../eur/chr()/hipstr.eur.vcf.gz'.format(chrom)
	hipstr_command = 'bcftools query -r ' + positions_arg + ' -f "%REF %ALT [%GT]\n" -s ' + sample + ' ' + hipstr_file_loc  #TODO the format notation
	hipstr_lines = sp.run(hipstr_command, shell = True, stdout = sp.PIPE, stderr = sp.PIPE).stdout.decode().read_lines()

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

		imputed_gp


for chrom in range(1,23):
	positions = find_overlapping_STRs(chrom)
	
	
	sample_list = get_sample_list() #TODO
	for sample in sample_list:
		gather_data(
	#for each ppt, gather those loci
	#analyze the results
	break

#make the plots
#plots are histograms of predicted accuracy vs actual accuracy
#with error bars denoting what the expected value is
#one plot with just the mostly likely predicted values (genotypes)
#one plot with the most likely predicted phased genotypes
#one plot with all predicted value
#one plot with only the non predicted values
#stratify some plots by the distribution of alleles in the overall population

