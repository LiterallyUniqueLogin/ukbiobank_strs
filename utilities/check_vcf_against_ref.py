import python_utils
import argparse
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument('--ref', required=True)
parser.add_argument('--vcf', required=True, nargs = '+')
args = parser.parse_args()

print("Loading reference", end='\r')
refs = python_utils.load_reference(args.ref)
bad_ref_variants = set()
for vcf in args.vcf:
	print(f"Working on vcf {vcf}", end = '\r')
	command = f'''source ~/.bashrc && conda activate bcftools && \
	bcftools query -f '%CHROM %POS %REF\n' {vcf}'''
	proc = sp.Popen(command, shell = True, stdout = sp.PIPE, stderr = sp.PIPE, universal_newlines = True)
	output = proc.stdout

	no_output = True
	for line_num, line in enumerate(output):
		no_output = False
		if line_num % 1000 == 0:
			print(f'Working on variant {line_num}', end ='\r')
		chr, pos, ref = line.split()
		pos = int(pos)
		recognized_chr = False
		for i in range(1,23):
			if chr == str(i) or chr == 'chr' + str(i):
				if not refs[i][(pos-1):(pos + len(ref) - 1)] == ref:
					bad_ref_variants.add((vcf, chr, pos))
				recognized_chr = True
				break
		if not recognized_chr:
			print(f"Couldn't recognize chr {chr} at pos {pos}")
			exit(-1)
	if no_output:
		print("Didn't read any lines from the vcf. Printing out bcftools stderr:")
		for line in proc.stderr:
			print(line)
		exit(-1)

if len(bad_ref_variants) > 0:
	print("The following variants do not match the reference")
	print('\n'.join(str(x) for x in sorted(bad_ref_variants)))
	exit(1)

print("All variants match the reference")
exit(0)
