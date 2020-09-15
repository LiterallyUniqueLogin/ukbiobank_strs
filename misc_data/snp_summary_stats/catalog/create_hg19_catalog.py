import csv

import numpy as np

hg19_bed = np.genfromtxt("catalog_hg19.bed", dtype=object)
#  Assert all positions are still 1bp wide
assert np.all(hg19_bed[:, 2].astype(int) - hg19_bed[:, 1].astype(int) == 1)

# convert 'chr9' to int(9)
hg19_chrom = np.char.rpartition(hg19_bed[:, 0].astype('U10'), 'chr')
assert not np.any(hg19_chrom[:, 1] == '')
hg19_chrom = hg19_chrom[:, 2].astype(int)
assert np.all(np.logical_and(1 <= hg19_chrom, hg19_chrom <= 22))

hg19_pos = hg19_bed[:, 2].astype(int)
row_nums = hg19_bed[:, 3].astype(int)


with open("catalog_hg38.tsv") as reader_file, \
        open("catalog_hg19.tsv", 'w') as writer_file:
    reader = csv.reader(reader_file, delimiter='\t')
    writer = csv.writer(writer_file, delimiter='\t')
    first = True
    for row_num, row in enumerate(reader):
        if first:
            writer.writerow(row)
            first = False
        print(f'row_num: {row_num}', end='\r')
        if row_num in row_nums:
            row[11] = int(hg19_chrom[row_nums == row_num])
            row[12] = int(hg19_pos[row_nums == row_num])
            writer.writerow(row)

