import glob
import os
import sys
import time
import merge_all_functions as merge

if len(sys.argv) != 2:
        print("Expecting exactly 1 argument, the chromosome number")
        exit(-1)
chr = sys.argv[1]

#set n_samples for all files and then last file
samples_per_file = 1000
samples_last_file = 409
n_files = len(glob.glob("{}/str_imputed/hap_no_preqc/vcf_batches/chr{}_samples_*.vcf" .format(\
	os.environ['UKB'], chr)))
n_files_minus_one = n_files - 1 #precompute this

try:
	#set up the overall script and open all the files
	startup_start = time.time()
	input_vcfs = []
	vcf_names = []
	output = None
	sample_block_len = len('\t0|0:0:0:0')
	for i in range(0, n_files):
		file_name = glob.glob("{}/str_imputed/hap_no_preqc/vcf_batches/chr{}_samples_{}_*.vcf" .format(\
			os.environ['UKB'], chr, i*samples_per_file + 1))[0]
		vcf_names.append(file_name)
		input_vcfs.append(open(file_name, 'br'))
	first_vcf = input_vcfs[0]

	output_loc = "{}/str_imputed/hap_no_preqc/vcfs/chr{}.vcf".format(os.environ['UKB'], chr)
	preexisting_output = os.path.exists(output_loc)

	if not preexisting_output:
		print("Starting the new merge")
		output = open(output_loc, 'wb')
		merge.validate_and_write_metadata(input_vcfs, output)
	else:
		print("Found a file already in progress, continuing from there")
		output = open(output_loc, 'r+b')
		variant_pos = merge.truncate_last_defined_variant(output)
		for i, vcf in enumerate(input_vcfs):
			merge.seek_to_variant(vcf, variant_pos, vcf_names[i])
			print('Done setting up file {}/{}   '.format(i, len(input_vcfs)), end='\r')
	#either way, all of the input_vcfs will be pointing to the next
	#character after the newline in front of the next variant to write
	#out to the output

	print("startup time: {:.2}s".format(time.time() - startup_start))
	
	#Write each variant to the vcf
	#assume all data lines are properly formatted
	all_variants_start = time.time()
	n_variants = 0 #number of variants already emitted

	while True: #iterate over all variants
		if n_variants % 10 == 0:
			variant_start = time.time()
		byte = first_vcf.read(1)
		fixed_field_bytes = 0
		if byte == b"": #check if we're done
			print("Checking if we're done")
			for i, vcf in enumerate(input_vcfs[1:]):
				i += 1
				byte = vcf.read(1)
				if byte != b"":
					print("First file {} is done but file {} isn't".format(vcf_names[0], vcf_names[i]))
					exit(-1)
			print("Done")
			exit(0)
		tab_count = 0
		info = bytearray()
		#read the fixed fields from the first file
		alt_alleles = 1
		while tab_count < 9:
			fixed_field_bytes += 1
			if byte == b'\t':
				tab_count += 1
				if tab_count == 8:
					if b'IMP' in info:
						output.write(b'\tIMP')
					else:
						output.write(b'\t.')

			if tab_count == 4 and byte == b',':
				alt_alleles += 1

			if tab_count != 7:
				output.write(byte)
			elif tab_count == 7:
				info += byte
				
			byte = first_vcf.read(1)
			if byte == b"":
				print("Didn't expect to be done! Finished before samples started In vcf {}".format(vcf_names[0]))
				exit(-1)
		output.write(byte)
		min_len = (4 + 6 * alt_alleles) * samples_per_file - 2
		#emit the sample data from the first file
		output.write(first_vcf.read(min_len))
		byte = first_vcf.read(1)
		while byte != b'\n':
			if byte == b'':
				print("Didn't expect to be done! Finished after samples started but before newline in vcf {}".format(vcf_names[0]))
				exit(-1)
			output.write(byte)
			byte = first_vcf.read(1)
		#emit the sample data from all other files (assume the fixed fields are correct)
		i = 0
		for vcf in input_vcfs[1:]:
			i += 1
			output.write(b'\t')
			#skip fixed fields. Assume that all fixed fields are exactly the same length
			#(seems to be true, all info fields are fixed precision)
			for _ in range(fixed_field_bytes):
				vcf.read(1)

			#output all sample data
			if i != n_files_minus_one:
				output.write(vcf.read(min_len))
			else:
				output.write(vcf.read((4 + 6 * alt_alleles) * samples_last_file - 2))
			byte = vcf.read(1)
			while byte != b'\n':
				if byte == b"":
					print("Didn't expect to be done! Ran into EOF before newline in vcf {}".format(vcf_names[i]))
					exit(-1)
				output.write(byte)
				byte = vcf.read(1)
				
		output.write(b'\n')

		now = time.time()
		n_variants += 1
		if n_variants % 1000 == 0:
			print('Variants complete(d this run): {}, Total time: {:.2}s, time/variant: {:.2}s, last 10 variants time: {:.2}s                   '.format(
				n_variants, now - all_variants_start, (now - all_variants_start)/n_variants, now - variant_start), end = '\r')
finally:
	for file in input_vcfs:
		if file is not None:
			file.close()
	if output is not None:
		output.close()
					
