import os

def produceListFromSampleQC(listName, keep, column):
	if keep:
		keepOrRemove = "keep"
	else:
		keepOrRemove = "remove"
	with open(os.environ["UKB"] + "/non_genetic_data/EGA/ukb_sqc_v2.txt") as csv, \
		open(os.environ["UKB"] + "/original/bgen_original/hap/ukb46122_cal_chr1_v2_s488282.fam") as fam, \
		open(os.environ["UKB"] + "/post_imputation_qc/sample_filtering/hap/" + keepOrRemove + "/" + listName + ".sample", 'w') as keep:
		with open(os.environ["UKB"] + "/original/bgen_original/hap/ukb46122_hap_chr1_v2_s487314.sample") as sample: 
			keep.write(sample.readline())

		#Assumes already validated taht csv and fam are the same length
		for csvline in csv:
			fam_tokens = fam.readline().split()
			if bool(int(csvline.split()[column])):
				keep.write("{} {} {} {}\n".format(fam_tokens[0], fam_tokens[1], fam_tokens[2], fam_tokens[4]))

		
