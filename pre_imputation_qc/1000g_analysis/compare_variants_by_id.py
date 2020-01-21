import os
import subprocess as sp

nid_overlap = 0
nid_no_pos_overlap = 0
ukb = os.environ['UKB']
datasets = os.environ['DATASETS']

for chr in range(1, 23):
	print(f"Working on chr{chr}, total so far {nid_overlap}", end='\r')

	ref_panel_ids = set()
	ref_panel_locs = {} #tuples (start, end) inclusive
	for line in sp.Popen(f"grep -v '^#' {ukb}/snpstr/vcf_1_sample/chr{chr}.vcf | cut -f 2-4", \
		shell = True, stdout=sp.PIPE, universal_newlines=True).stdout:
		pos, id, ref = line.split()
		ref_panel_ids.add(id)
		ref_panel_locs[id] = (int(pos), int(pos) + len(ref) - 1)

	hipstr_ids = set()
	hipstr_locs = {} #tuples (start, end) inclusive
	for line in sp.Popen(f"source ~/.bashrc && \
			conda activate bcftools && \
			bcftools query -f '%POS %ID %REF\n' \
			{datasets}/1000Genomes/hipstr_calls_sample_subset/eur/chr{chr}/hipstr.chr{chr}.eur.vcf.gz", \
		shell = True, stdout=sp.PIPE, universal_newlines=True).stdout:
		pos, id, ref = line.split()
		hipstr_ids.add(id)
		hipstr_locs[id] = (int(pos), int(pos) + len(ref) - 1)

	overlap = hipstr_ids.intersection(ref_panel_ids)
	with open(f"{ukb}/pre_imputation_qc/1000g_analysis/output/variant_id_overlaps/chr{chr}/overlap.ids", 'w') as overlap_file, \
		open(f"{ukb}/pre_imputation_qc/1000g_analysis/output/variant_id_overlaps/chr{chr}/id_overlap_but_not_seq_overlap.ids", 'w') as not_seq_file:
		for ID in sorted(overlap):
			overlap_file.write(ID)
			overlap_file.write("\n")
			lesser = min(ref_panel_locs[ID], hipstr_locs[ID])
			greater = max(ref_panel_locs[ID], hipstr_locs[ID])
			if lesser[1] < greater[0]:
				not_seq_file.write(ID)
				not_seq_file.write('\n')
				nid_no_pos_overlap += 1
	nid_overlap += len(overlap)

print(f"""The hipstr calls overlap the reference panel at {nid_overlap} loci,
	as determined solely by locus id. Of those, {nid_no_pos_overlap} don't overlap at all
	in reference sequence location.""")
