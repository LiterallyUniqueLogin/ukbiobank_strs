import argparse
import os
import os.path
import re
import subprocess as sp
import sys

import matplotlib.pyplot as plt
import numpy as np

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('summary_file_name')
args = parser.parse_args()

if not os.path.exists(args.summary_file_name):
    print(f'File {args.summary_file_name} does not exist')
    sys.exit(1)

str_lens = {}

with open(f'{ukb}/side_analyses/exome_strs/snpstr_exome_strs_38.bed') as exome_strs:
    for line in exome_strs:
        split = line.split()
        str_lens[split[5]] = int(split[2]) - int(split[1])

max_called_lens = {}

for chrom in range(1, 23):
    with open(f'{ukb}/side_analyses/exome_strs/snpstr_panel_exome_calls_hg19_chr{chrom}.txt') as calls:
        for str_locus in calls:
            str_locus = str_locus.strip()
            if str_locus == '':
                continue
            idee, ref, alt, gts = str_locus.split()
            allele_lens = np.array(
                [len(allele) - len(ref) + str_lens[idee] for allele in [ref] + alt.split(',')]
            )
            alleles = np.array([int(x) for x in re.split(r':|\||/', gts) if x])
            max_called_lens[idee] = max(np.max(allele_lens[alleles]),
                                        len(ref))

assert len(max_called_lens) == 630

min_max_len = np.min(np.array([x for x in max_called_lens.values()]))
max_max_len = np.max(np.array([x for x in max_called_lens.values()]))

failed_strs = {}
failed_str_batches = set()
with open(args.summary_file_name) as summary:
    for line in summary:
        parts = line.split('_batch_')
        if len(parts) < 2:
            continue
        if line in failed_str_batches:
            continue
        failed_str_batches.add(line)
        STR = parts[0]
        if STR not in failed_strs:
            failed_strs[STR] = 1
        else:
            failed_strs[STR] += 1

filtered_ids = set()
with open('filtered_ids.txt') as filtered_ids_file:
    for line in filtered_ids_file:
        line = line.strip()
        if '#' in line or line == '':
            continue
        filtered_ids.add(line)

successes = []
fails = []
for STR, max_len in max_called_lens.items():
    if STR in filtered_ids:
        fails += [max_len] * 10
        continue
    if STR in failed_strs:
        n_fails = failed_strs[STR]
        assert n_fails <= 10
        fails += [max_len] * n_fails
        successes += [max_len] * (10 - n_fails)
    else:
        successes += [max_len]*10

print((len(successes) + len(fails))/10)

plt.hist([successes, fails], range(min_max_len-1, max_max_len+2),
         stacked=True, histtype='bar', label=['successes', 'fails'],
         color=['green', 'red'])
plt.legend()
plt.ylabel('# loci')
plt.xlabel('max called len at locus in SNPSTR panel')
plt.savefig(f'{args.summary_file_name}.distribution.png')
