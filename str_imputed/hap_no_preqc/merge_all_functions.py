import datetime

class MergeException(Exception):
	pass

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
			raise MergeException()

		for line1, line2 in zip(metadata_lines, metadata):
			if not ('filedate' in line1 and 'filedate' in line2) and line1 != line2:
				print('Found file {} with different metadata than file {}'.format(vcf_names[i], vcf_names[0]))
				print('line1', line1)
				print()
				print('line2', line2)
				raise MergeException()
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
			raise MergeException()

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
	neg_file_size = - file.seek(0, 2)
	buffer_len = 2**16
	file.seek(-buffer_len, 2)
	buffer = file.read(buffer_len)
	buffer_idx = 0
	steps_back = 0
	jumps_back = 1
	#Find the last newline
	while True:
		if buffer_idx == buffer_len:
			jumps_back += 1
			buffer_idx = 1
			next_pos = -jumps_back*buffer_len
			if next_pos < neg_file_size:
				print("Found no variants in the preexisting output file!")
				raise MergeException()
			file.seek(next_pos, 2)
			buffer = file.read(buffer_len)
		else:
			buffer_idx += 1
		byte = buffer[(buffer_len-buffer_idx):(buffer_len-buffer_idx + 1)]
		steps_back += 1
		if byte == b'\n':
			file.seek(-buffer_idx-(jumps_back-1)*buffer_len, 2)
			#Check if this is a completely defined index or not
			tab_count = 0
			contents = bytearray()
			byte = file.read(1)
			while byte != b'':
				if byte == b'\t':
					tab_count += 1
				elif tab_count == 1:
					contents += byte
				if tab_count == 2:
					#terminating condition
					file.seek(1-steps_back, 2)
					file.truncate()
					return contents
				byte = file.read(1)
			#hit end of file before hitting the second tab, need to go to the previous newline
			file.seek(-steps_back-1, 2)
			#go back to the overall search loop

#returns a tuple with the variant position in bytes (e.g. b'92243')
#and the number of bytes read to get to the current location
#or (None, None) if end of file is reached
#don't read past bound while doing this
def read_till_next_variant_pos(vcf, start_idx, bound):
	bound += 50 #make sure we don't set the bound in the middle of the variant name at the start of the line
	steps_left = bound - start_idx 
	steps_taken = 0
	byte = b'a' #dummy value
	newline_encountered = False
	tabs_encountered = 0
	pos = bytearray()
	while byte != b'' and steps_left >= 0:
		steps_left -= 1
		steps_taken += 1
		byte = vcf.read(1)
		if not newline_encountered:
			#keep reading till the beginning of the next line
			if byte == b'\n':
				newline_encountered = True
			continue
		if byte == b'#':
			#skip header or metadata lines
			newline_encountered = False
			continue
		if tabs_encountered == 0:
			#skip the chrom number column
			if byte == b'\t':
				tabs_encountered = 1
			continue
		if tabs_encountered == 1:
			#read off the position column
			if byte == b'\t':
				return (pos, steps_taken)
			else:
				pos += byte
	return (None, None) #hit the bound (which might be EOF) before finding the variant

#navigates to the beginning of the line denoted by the variant pos
#(e.g. b'482390')
#and returns the file position of the line's first byte
def seek_to_variant(vcf, variant_pos, vcf_name):
	#binary search
	min = 0
	max = vcf.seek(0, 2)
	while True:
		half = (min + max)//2
		vcf.seek(half)
		pos, steps_taken = read_till_next_variant_pos(vcf, half, max)
		if pos == None:
			print("Couldn't find variant with pos {} in file {}".format(variant_pos.decode('us-ascii'), vcf_name))
			raise MergeException()
		elif pos == variant_pos:
			current_idx = half + steps_taken
			break
		elif len(pos) < len(variant_pos) or (len(pos) == len(variant_pos) and pos < variant_pos):
			min = half + steps_taken
			continue
		else:
			max = half
			continue
	#seek backwards till we get to just after the newline
	byte = vcf.read(1)
	while byte != b'\n':
		current_idx -= 1
		vcf.seek(current_idx)
		byte = vcf.read(1)
	return current_idx + 1
	
#functions for debugging
def print_till_EOF(file):
	all_bytes = bytearray()
	byte = file.read(1)
	while byte != b'':
		all_bytes += byte
		byte = file.read(1)
	print(all_bytes.decode('us-ascii'))

def print_till_newline(file):
	all_bytes = bytearray()
	byte = file.read(1)
	while byte != b'' and byte != b'\n':
		all_bytes += byte
		byte = file.read(1)
	print(all_bytes.decode('us-ascii'))

