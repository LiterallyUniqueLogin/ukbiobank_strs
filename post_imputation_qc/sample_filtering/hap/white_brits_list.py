import os

with open(os.environ["UKB"] + "/non_genetic_data/EGA/ukb_sqc_v2.txt") as csv, \
	open(os.environ["UKB"] + "/original/bgen_original/hap/ukb46122_cal_chr1_v2_s488282.fam") as fam, \
	open(os.environ["UKB"] + "/post_imputation_qc/sample_filtering/hap/keep/white_british.sample", 'w') as keep:

	with open(os.environ["UKB"] + "/original/bgen_original/hap/ukb46122_hap_chr1_v2_s487314.sample") as sample: 
		keep.write(sample.readline())

	#Assumes already validated taht csv and fam are the same length
	for csvline in csv:
		fam_tokens = fam.readline().split()
		if bool(int(csvline.split()[23])):
			keep.write("{} {} {} {}\n".format(fam_tokens[0], fam_tokens[1], fam_tokens[2], fam_tokens[4]))

	
