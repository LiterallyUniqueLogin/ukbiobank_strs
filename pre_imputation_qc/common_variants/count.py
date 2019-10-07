import subprocess as sp

with open('counts.csv', 'w') as outfile:
	outfile.write('chr, #SNPSTR Ref Panel Varaints, #UKB Directly Genotyped Variants, #Overlap, %UKB lost\n')

	#run with conda activate bcftools
	for i in range(1,23):
		count_hap_cmd = 'wc -l $UKB/original/vcf_1_sample/hap/chr{}.vcf'.format(i)
		done = sp.run(count_hap_cmd, shell = True, capture_output = True).stdout.decode()
		hap_count = int(done.split()[0]) - 6
		count_snpstr_cmd = 'wc -l $UKB/snpstr/vcf_1_sample/chr{}.vcf'.format(i)
		done = sp.run(count_snpstr_cmd, shell = True, capture_output = True).stdout.decode()
		snpstr_count = int(done.split()[0]) - 11
		isec_cmd = './intersect.sh {} | wc -l'.format(i)
		overlap = int(sp.run(isec_cmd, shell = True, capture_output = True).stdout.decode())
		outfile.write('{}, {}, {}, {}, {:.2%}\n'.format(i, snpstr_count,
			hap_count, overlap, 1 - overlap/hap_count))



