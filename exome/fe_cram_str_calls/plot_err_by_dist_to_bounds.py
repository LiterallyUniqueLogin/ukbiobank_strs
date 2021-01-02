import os

import matplotlib.pyplot as plt
import numpy as np

ukb = os.environ['UKB']

filtered_ids = set()
with open('filtered_ids.txt') as filtered_ids_file:
    for line in filtered_ids_file:
        line = line.strip()
        if '#' in line or line == '':
            continue
        filtered_ids.add(line)

str_pos = {}

with open(f'{ukb}/side_analyses/exome_strs/snpstr_exome_strs_38.bed') as exome_strs:
    for line in exome_strs:
        split = line.split()
        str_pos[split[5]] = (int(split[0][3:]), int(split[1]), int(split[2]))

exome_chr_start_locs = {}
exome_chr_stop_locs = {}

for chrom in range(1, 23):
    exome_chr_start_locs[chrom] = []
    exome_chr_stop_locs[chrom] = []

str_chroms = set([str(x) for x in range(1,23)])

with open(f'{ukb}/side_analyses/exome_strs/intermediate_files/exome_38.bed') as exome:
    for line in exome:
        split = line.split()
        if split[0] not in str_chroms:
            continue
        exome_chr_start_locs[int(split[0])].append(int(split[1]))
        exome_chr_stop_locs[int(split[0])].append(int(split[2]))

for chrom in range(1, 23):
    exome_chr_start_locs[chrom] = np.array(exome_chr_start_locs[chrom])
    exome_chr_stop_locs[chrom] = np.array(exome_chr_stop_locs[chrom])

successes = []
fails = []

for STR in str_pos:
    chrom, start, stop = str_pos[STR]
    start_locs = exome_chr_start_locs[chrom]
    start_dist = np.min((start - start_locs)[start_locs <= start])
    stop_locs = exome_chr_stop_locs[chrom]
    stop_dist = np.min((stop_locs - stop)[stop <= stop_locs])
    dist = np.log(min(start_dist, stop_dist) + 1)
    if STR in filtered_ids:
        fails.append(dist)
    else:
        successes.append(dist)

min_dist = min(min(successes), min(fails))
max_dist = max(max(successes), max(fails))

plt.hist([successes, fails],
         stacked=True, histtype='bar', label=['successes', 'fails'],
         color=['green', 'red'])
plt.legend()
plt.ylabel('# loci')
plt.xlabel('log distance to exome capture boundary')
plt.savefig('dist_distribution.png')
