import os
import subprocess as sp

nid_overlap = 0
ukb = os.environ['UKB']
datasets = os.environ['DATASETS']

for chr in range(1, 23):
	print(f"Working on chr{chr}, total so far {nid_overlap}", end='\r')

	ref_panel_ids = set()
	for line in sp.run(f"grep -v '^#' {ukb}/snpstr/vcf_1_sample/chr{chr}.vcf | cut -f 3", \
		shell = True, stdout=sp.PIPE, universal_newlines=True).stdout.split('\n'):
		ref_panel_ids.add(line.strip())

	hipstr_ids = set()
	for line in sp.run(f"zcat {datasets}/1000Genomes/hipstr_calls_sample_subset/eur/chr{chr}/hipstr.chr{chr}.eur.vcf.gz | grep -v '^#' | cut -f 3", \
		shell = True, stdout=sp.PIPE, universal_newlines=True).stdout.split('\n'):
		hipstr_ids.add(line.strip())

	nid_overlap += len(hipstr_ids.intersection(ref_panel_ids))

print(f"""The hipstr calls overlap the reference panel at {nid_overlap} loci,
	as determined solely by locus id""")
