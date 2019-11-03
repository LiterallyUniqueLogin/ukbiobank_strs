import glob
import os
import sys
import datetime
import time

if len(sys.argv) != 2:
        print("Expecting exactly 1 argument, the chromosome number")
        exit(-1)
chr = sys.argv[1]

def validate_and_write_metadata(input_vcfs, output):
	metadata_lines = bytearray()
	metadata_to_write = ""

	first_vcf = input_vcfs[0]
	byte = first_vcf.read(1)
	just_saw_newline = True
	num_hashes_after_newline = 0
	while True: #read the metadata till the #CHROM line
		metadata_lines += byte
		
		if just_saw_newline:
			if byte == b'#':
				num_hashes_after_newline += 1
			elif num_hashes_after_newline < 2:
				break
			else:
				num_hashes_after_newline = 0
				just_saw_newline = False
		else:
			just_saw_newline = byte == b'\n'
		byte = first_vcf.read(1)

	byte = first_vcf.read(1)
	tab_count = 0
	while True: #Read till the first sample ID
		metadata_lines += byte
		if byte == b'\t':
			tab_count += 1
			if tab_count == 9:
				break
		byte = first_vcf.read(1)

	#Confirm all other headers are the same except for maybe the date line
	metadata_len = len(metadata_lines)
	metadata_lines = metadata_lines.decode('us-ascii').split('\n')

	for i, file in enumerate(input_vcfs[1:]):
		i += 1
		metadata = file.read(metadata_len)
		metadata = metadata.decode('us-ascii').split('\n')
		if len(metadata) != len(metadata_lines):
			print('File {} has obviously different metadata than file {}'.format(vcf_names[i], vcf_names[0]))
			exit(-1)

		for line1, line2 in zip(metadata_lines, metadata):
			if not ('filedate' in line1 and 'filedate' in line2) and line1 != line2:
				print('Found file {} with different metadata than file {}'.format(vcf_names[i], vcf_names[0]))
				print('line1', line1)
				print()
				print('line2', line2)
				exit(-1)
	print('All files have the same metadata')

	#Write the correct metadata to the output file
	for line in metadata_lines:
		bin_line = (line + "\n").encode('us-ascii') #binary encoding of the line
		#If modifying these, also modify the fixed fields section below
		if '##fileformat' in line:
			output.write(bin_line)
		elif '##filedate' in line:
			output.write('##filedate={}\n'.format(datetime.datetime.now().strftime('%Y%m%d')).encode('us-ascii'))
		elif '##source' in line:
			output.write('##source=merge_all.py\n'.encode('us-ascii'))
			output.write(bin_line)
		elif '##INFO=<ID=AF' in line:
			continue
		elif '##INFO=<ID=DR2' in line:
			continue
		elif '##INFO=<ID=IMP' in line:
			output.write(bin_line)
		elif '##FORMAT=' in line:
			output.write(bin_line)
		elif '#CHROM' in line:
			output.write(line.encode('us-ascii'))
		else:
			print('Unrecognized metadata line')
			exit(-1)

	#Write all the sample ids to the header line
	first = True
	for i, vcf in enumerate(input_vcfs):
		if not first:
			output.write(b'\t')
		first = False
		
		byte = vcf.read(1)
		while byte != b'\n':
			output.write(byte)
			byte = vcf.read(1)
	output.write(b'\n')

#seeks to the point in the file right after the \n before the last
#chrNum\tpos\t
#and trucates the file.
#Returns pos
def truncate_last_defined_variant(file):
        buffer_len = 2**16
        file.seek(-buffer_len, 2)
        buffer = file.read(buffer_len)
        buffer_idx = 0
	#Find the next newline
	while True:
                if buffer_idx == buffer_len:
                        buffer_idx = 1
                        file.seek(-2*buffer_len, 1)
                        buffer = file.read(buffer_len)
                else:
                        buffer_idx += 1
                byte = buffer[(buffer_len-buffer_idx):(buffer_len-buffer_idx + 1)]
                if byte == b'\n':
			#Check if this is a completely defined index or not
                	file.seek(1-buffer_idx, 1) #this puts you right after the newline
			line_start = file.tell()
			tab_count = 0
			contents = bytearray()
			byte = file.read(1)
			while byte != b'':
				contents += byte
				if byte == b'\t':
					tab_count += 1
				if tab_count == 2:
					#terminating condition
					return_pos = contents.decode('us-ascii').split()[1]
					file.seek(line_start)
					file.truncate()
					return return_pos
			#hit end of file before hitting the second tab, need to go to the previous newline
			file.seek(line_start - 1, 0)
			#go back to the overall search loop

def seek_to_variant(vcf, variant_pos):
	

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

	output_loc = "{}/str_imputed/hap_no_preqc/vcfs/chr{}.vcf".format(os.environ['UKB'], chr)
	preexisting_output = os.path.exists(output_loc)
	output = open(output_loc, 'r+b')

	if not preexsiting_output:
		validate_and_write_metadata(input_vcfs, output)
	else:
		variant_pos = truncate_last_defined_variant(output)
		for vcf in input_vcfs:
			seek_to_variant(vcf, variant_pos)
	#either way, all of the input_vcfs will be pointing to the next
	#character after the newline in front of the next variant to write
	#out to the output

	print("startup time: {:.2}s".format(time.time() - startup_start))
	
	#Write each variant to the vcf
	#assume all data lines are properly formatted
	all_variants_start = time.time()
	n_variants = 0 #number of variants already emitted

	#minimum length of the data portion of a line
	min_len = samples_per_file*sample_block_len - 2
	min_len_last_file = samples_last_file*sample_block_len - 2
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
		while tab_count < 9:
			fixed_field_bytes += 1
			if byte == b'\t':
				tab_count += 1
				if tab_count == 8:
					if b'IMP' in info:
						output.write(b'IMP')
					else:
						output.write(b'.')
				output.write(byte)
			elif tab_count < 7 or tab_count == 8:
				output.write(byte)
			elif tab_count == 7:
				info += byte
				
			byte = first_vcf.read(1)
			if byte == b"":
				print("Didn't expect to be done! Finished before samples started In vcf {}".format(vcf_names[0]))
				exit(-1)
		#emit the sample data from the first file
		output.write(byte)
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
			vcf.seek(fixed_field_bytes, 1)

			#output all sample data
			if i != n_files_minus_one:
				output.write(vcf.read(min_len))
			else:
				output.write(vcf.read(min_len_last_file))
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
		if n_variants % 10 == 0:
			print('Variants complete: {}, Total time: {:.2}s, time/variant: {:.2}s, last 10 variants time: {:.2}s                   '.format(
				n_variants, now - all_variants_start, (now - all_variants_start)/n_variants, now - variant_start), end = '\r')

finally:
	for file in input_vcfs:
		if file is not None:
			file.close()
	if output is not None:
		output.close()
					
