import bgen_reader

print("Reading fasta", flush=True)
with open('/projects/ps-gymreklab/resources/dbase/human/hg19/hg19.fa') as fasta:
    genome_string = fasta.read()
print("Munging fasta", flush=True)
genome_string = genome_string.replace('\n', '')
genome_string = genome_string.upper()
genome_list = genome_string.split('>')
chroms = {}
for chrom_seq in genome_list:
    if len(chrom_seq) < 4:
        continue
    if '_' in chrom_seq[:100]:
        continue
    print(chrom_seq[0:10])
    if chrom_seq[4] in '0123456789':
        chroms[int(chrom_seq[3:5])] = chrom_seq[5:]
    else:
        if chrom_seq[3] not in '0123456789':
            continue
        chroms[int(chrom_seq[3])] = chrom_seq[4:]

assert len(chroms) == 22
for chrom in range(1, 23):
    assert chrom in chroms

for chrom in range(1, 23):
    print(f"Working on chrom {chrom}", flush=True)
    with bgen_reader.open_bgen(f'ukb_imp_chr{chrom}_v3.bgen',
                               allow_complex=True,
                               verbose=False) as bgen:
        for nvar, (pos, alleles) in enumerate(zip(bgen.positions,
                                                  bgen.allele_ids)):
            if nvar % 10000 == 0:
                print(f"Working on variant {nvar}", end='\r', flush=True)
            presumed_ref = alleles.split(',')[0]
            actual_ref = chroms[chrom][(pos-1):(pos+len(presumed_ref)-1)]
            if actual_ref != presumed_ref:
                print(f"{chrom}:{pos} expected {actual_ref} got {alleles}",
                      flush=True)


