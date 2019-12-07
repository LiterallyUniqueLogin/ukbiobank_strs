import os
import subprocess as sp

nid_overlap = 0
ukb = os.environ['UKB']
datasets = os.environ['DATASETS']

for chr in range(1, 23):
	print(f"Working on chr{chr}, total so far {nid_overlap}", end='\r')

	ref_panel_ids = set()
	for line in sp.Popen(f"grep -v '^#' {ukb}/snpstr/vcf_1_sample/chr{chr}.vcf | cut -f 3", \
		shell = True, stdout=sp.PIPE, universal_newlines=True).stdout:
		ref_panel_ids.add(line.strip())

	hipstr_ids = set()
	for line in sp.Popen(f"source ~/.bashrc && \
			conda activate bcftools && \
			bcftools query -f '%ID\n' \
			{datasets}/1000Genomes/hipstr_calls_sample_subset/eur/chr{chr}/hipstr.chr{chr}.eur.vcf.gz", \
		shell = True, stdout=sp.PIPE, universal_newlines=True).stdout:
		hipstr_ids.add(line.strip())

	overlap = hipstr_ids.intersection(ref_panel_ids)
	with open(f"{ukb}/pre_imputation_qc/ref_panel_size/output/variant_id_overlaps/chr{chr}/overlap.ids", 'w') as overlap_file:
		for ID in sorted(overlap):
			overlap_file.write(ID)
			overlap_file.write("\n")
	nid_overlap += len(overlap)

print(f"""The hipstr calls overlap the reference panel at {nid_overlap} loci,
	as determined solely by locus id""")
