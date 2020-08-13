import os

ukb = os.environ['UKB']
basedir = f"{ukb}/side_analyses/entropy"

for chrom in range(1, 23):
    with open(f"{basedir}/full_genome/chr{chrom}.tab") as full_genome, \
            open(f"{basedir}/exome/chr{chrom}.tab") as exome, \
            open(f"{basedir}/combined/chr{chrom}.tab", 'w') as combined:
        genome_iter = iter(full_genome)
        exome_iter = iter(exome)
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
