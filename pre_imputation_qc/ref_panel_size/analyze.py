import subprocess as sp

sample = 'HG00351'
chr = 1

command = """
source ~/.bashrc && \
conda activate bcftools && \
bcftools query -s {} -f '%POS %ID %REF %ALT [%GT]' \
/projects/ps-gymreklab/resources/datasets/1000Genomes/hipstr_calls_sample_subset/eur/chr{}/hipstr.chr{}.eur.vcf.gz'
""".format(sample, chr, chr)

print(command)
