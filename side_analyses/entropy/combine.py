import argparse
import os

"""
combines the files in full_genome/ and exome/ into combined/
with two cols entropy-whole-genome and entropy-exome
most of the entries in the later are nan (for loci not in the exome)
"""

ukb = os.environ['UKB']
basedir = f"{ukb}/side_analyses/entropy"

parser = argparse.ArgumentParser()
parser.add_argument('--het', action='store_true')
args = parser.parse_args()
het = ""
if args.het:
    het = "_het"

for chrom in range(1, 23):
    with open(f"{basedir}/full_genome{het}/chr{chrom}.tab") as full_genome, \
            open(f"{basedir}/exome_filtered{het}/chr{chrom}.tab") as exome, \
            open(f"{basedir}/combined{het}/chr{chrom}.tab", 'w') as combined:
        genome_iter = iter(full_genome)
        exome_iter = iter(exome)
        genome_line = next(genome_iter)
        exome_line = next(exome_iter)
        #header
        combined.write(genome_line[:-1])
        if not args.het:
            combined.write("-whole-genome\tentropy-exome\n")
        else:
            combined.write("-whole-genome\thet-exome\n")
        genome_line = next(genome_iter)
        exome_line = next(exome_iter)

        _, genome_start, genome_end, _ = genome_line.split()
        _, exome_start, exome_end, exome_entropy = exome_line.split()
        while True:
            if genome_start == exome_start and genome_end == exome_end:
                combined.write(genome_line[:-1])
                combined.write(f"\t{exome_entropy}\n")
                try:
                    genome_line = next(genome_iter)
                    _, genome_start, genome_end, _ = genome_line.split()
                except StopIteration:
                    break

                try:
                    exome_line = next(exome_iter)
                    _, exome_start, exome_end, exome_entropy = exome_line.split()
                except StopIteration:
                    pass
            else:
                combined.write(genome_line[:-1])
                combined.write(f"\tnan\n")
                try:
                    genome_line = next(genome_iter)
                    _, genome_start, genome_end, _ = genome_line.split()
                except StopIteration:
                    break
