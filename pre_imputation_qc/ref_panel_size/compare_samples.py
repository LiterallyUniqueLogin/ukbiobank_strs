import subprocess as sp
import os

#Figure out which samples are in the hipstr call set but 
#not in the reference panel

samples_loc = os.environ['UKB'] + \
	"/pre_imputation_qc/ref_panel_size/output/compare_samples/"

command = """
plink2 --pfile $UKB/snpstr/pfile/chr1 \
	--write-samples \
	--out """ + samples_loc + "ref"

sp.run(command, shell=True)

ref_samples = set()
hipstr_samples = set()
hipstr_samples_to_ancestry = {}

with open(samples_loc + "ref.id") as ref_samples_file:
	first = True
	for line in ref_samples_file:
		if first:
			first = False
			continue
		ref_samples.add(line.split()[0])

for ancestry in ['eas', 'afr', 'eur']:
	command = """
zcat /projects/ps-gymreklab/resources/datasets/1000Genomes/hipstr_calls_sample_subset/{}/chr1/hipstr.chr1.{}.vcf.gz | grep CHROM | head -n 1
""".format(ancestry, ancestry)
	results = sp.run(command, shell = True, stdout = sp.PIPE)
	stdout = results.stdout.decode()
	samples = stdout.split()[9:]
	for sample in samples:
		hipstr_samples.add(sample)
		hipstr_samples_to_ancestry[sample] = ancestry

ancestry_missing = {'eas': 0, 'afr': 0, 'eur': 0}
missing_samples = set()
extant_eur_samples = set()
for sample in hipstr_samples:
	if sample not in ref_samples:
		ancestry_missing[hipstr_samples_to_ancestry[sample]] += 1
		missing_samples.add(sample)
	elif hipstr_samples_to_ancestry[sample] == 'eur':
		extant_eur_samples.add(sample)

with open(samples_loc + "hipstr_not_ref.samples", 'w') as missing_samples_file:
	missing_samples_file.write("ID\n")
	for sample in missing_samples:
		missing_samples_file.write(sample)
		missing_samples_file.write("\n")

with open(samples_loc + "eur_hipstr_and_ref.samples", 'w') as extant_samples_file:
	extant_samples_file.write("ID\n")
	for sample in extant_eur_samples:
		extant_samples_file.write(sample)
		extant_samples_file.write("\n")



print("""Found the following number of samples (per ancestry) that are in the
hipstr calls but not in the reference panel""")
print(ancestry_missing)

#Confirm that there are no samples in our reference panel
#that are not in the provided sample list 

sample_list_samples = set()
with open(os.environ['UKB'] + "/snpstr/igsr_samples.tsv") as sample_file:
	first = True
	for line in sample_file:
		if first:
			first = False
			continue
		sample = line.split()[0]
		sample_list_samples.add(sample)

for sample in ref_samples:
	if sample not in sample_list_samples:
		print("Found sample {} in the ref panel not in the sample list provided by 1000G".format(sample))
		exit(-1)

print("All samples in the reference panel are accounted for in the igsr sample list")
