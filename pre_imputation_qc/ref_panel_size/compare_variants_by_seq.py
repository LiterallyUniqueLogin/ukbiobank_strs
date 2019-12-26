import os
import subprocess as sp

ukb = os.environ['UKB']
datasets = os.environ['DATASETS']

def load_reference():
	with open('/projects/ps-gymreklab/resources/dbase/human/hg19/hg19.fa') as ref_file:
		lines = ref_file.readlines()
	contig_breaks = []
	contig_names = []
	for line_num, line in enumerate(lines):
		if line.startswith('>'):
			contig_breaks.append(line_num)
			contig_names.append(line[1:-1])
	refs = {}
	for contig_num, name in enumerate(contig_names):
		for chr in range(1,23):
			if name == 'chr{}'.format(chr):
				start_line_num = contig_breaks[contig_num] + 1
				end_line_num = contig_breaks[contig_num + 1]
				refs[chr] = (''.join(lines[start_line_num:end_line_num])).replace('\n', '')
				break
	return refs

#returns a tuple
#(hipstr_allele_lens, hipstr_allele_map, ref_allele_lens)
# each of this is two nested dics (e.g. hisptr_allele_lens[chr][ID])
# where chr is a string (e.g. '19')
# and ID is the ID of the variant (e.g. 'STR_259')
# hipstr_allele_lens : a list containing the length of the alleles relative to the reference (so (0, len(alt1) - len(ref), ... , len(altn) - len(ref)))
# hipstr_allele_map : mapping alleles numbers from the hipstr panel to their corresponding number in the reference panel (or to -1 if they don't occur in the reference panel)
# ref_allele_lens : like hipstr_allele_lens

def generate_overlaps
	hipstr_allele_lens = {}
	hipstr_allele_map = {}
	ref_allele_lens = {}

	noverlaps = 0

	for chr in range(1, 23):
		chr = str(chr)
		print(f"Working on chr{chr}", end='\r')
		hipstr_ids = {}
		ref_panel_ids = {}

		hipstr_allele_lens[chr] = {}
		hipstr_allele_map[chr] = {}
		ref_allele_lens[chr] = {}

		for line in sp.Popen(f"source ~/.bashrc && \
				conda activate bcftools && \
				bcftools query -f '%ID %REF,%ALT\n' \
				{ukb}/snpstr/1kg.snp.str.chr{chr}.vcf.gz", \
			shell = True, stdout=sp.PIPE, universal_newlines=True).stdout:
			line = line.strip().split()
			ref_panel_ids[chr][line[0]] = line[1].split(',')

		for line in sp.Popen(f"source ~/.bashrc && \
				conda activate bcftools && \
				bcftools query -f '%ID %REF,%ALT\n' \
				{datasets}/1000Genomes/hipstr_calls_sample_subset/eur/chr{chr}/hipstr.chr{chr}.eur.vcf.gz", \
			shell = True, stdout=sp.PIPE, universal_newlines=True).stdout:
			line = line.strip().split()
			hipstr_ids[chr][line[0]] = line[1].split(',')

		overlap = set(hipstr_ids.keys()).intersection(set(ref_panel_ids.keys()))
		noverlaps += overlap
		with open(f"{ukb}/pre_imputation_qc/ref_panel_size/output/variant_id_overlaps/chr{chr}/overlap.ids", 'w') as overlap_file:
			for ID in sorted(overlap):
				overlap_file.write(ID)
				overlap_file.write("\n")
		nid_overlap += len(overlap)

	print(f"""The hipstr calls overlap the reference panel at {nid_overlap} loci,
		as determined solely by locus id""")
