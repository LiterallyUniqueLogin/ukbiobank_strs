import merge_all_functions as merge
import os
import os.path
import shutil

test_dir = os.environ['UKB'] + "/str_imputed/hap_no_preqc/tests"
test_vcfs = ['not_pretruncated.vcf', 
	'pretruncated.vcf', 
	'pretruncated2.vcf',
	'pretruncated3.vcf',
	'pretruncated4.vcf']

print("Testing seek to reasonable variants")
with open(test_dir + "/" + 'not_pretruncated.vcf', 'br') as test_file:
	for variant_pos in ['92774', '60123']:
		variant_pos = variant_pos.encode('us-ascii')
		file_pos = merge.seek_to_variant(test_file, variant_pos, 'temp.vcf')
		merge.print_till_newline(test_file)
		test_file.seek(file_pos) #confirm that the returned file_pos is directly after the newline
		merge.print_till_newline(test_file)
		test_file.seek(-1000, 1)

print("Testing seek to truncated variant")
variant_pos = '99138'.encode('us-ascii')
for file_name in test_vcfs:
	with open(test_dir + "/" + file_name, 'br') as test_file:
		try:
			file_pos = merge.seek_to_variant(test_file, variant_pos, file_name)
		except merge.MergeException:
			pass
		merge.print_till_newline(test_file)
		test_file.seek(file_pos) #confirm that the returned file_pos is directly after the newline
		merge.print_till_newline(test_file)
		test_file.seek(-1000, 1)

print("Testing truncation")
for file_name in test_vcfs:
	dest_loc = test_dir + "/test_output_" + file_name
	shutil.copyfile(test_dir + "/" + file_name, dest_loc)
	with open(dest_loc, 'br+') as output_file:
		variant = merge.truncate_last_defined_variant(output_file)
		print(file_name, variant)
		merge.print_till_EOF(output_file)
		output_file.seek(-1000, 1)

