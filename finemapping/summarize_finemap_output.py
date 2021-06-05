#!/usr/bin/env python3

import argparse
import os
import os.path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('phenotype')
parser.add_argument('--signals', nargs='+')
args = parser.parse_args()
phenotype = args.phenotype
signals = ['/'.join(signal.split('/')[:-1]) for signal in args.signals]

# ---- what do I want to calculate
# ranking of most significant STRs across all loci
# percent and raw count loci with an STR in the most probable configuration, and list of such
# precent and raw count loci where STR is most probable, and list of such
# percent and raw count STRs that reach a certain threshold (ROC curve)
# percent and raw count SNPs that reach a certain threshold (ROC curve)
# loci ranked by highest probability - call out loci with low prob
# loci with number of causal SNPs - loci with too many causal SNPs
# reproducability?

# number of causal variants per locus

max_causal = 20

total_causal_count = np.zeros((max_causal))

max_causal_loci = set()

n_signals_with_strs = 0
for signal in signals:
    if os.path.exists(f'{signal}/no_strs'):
        continue
    n_signals_with_strs += 1
    found_cred_file = False
    for n in range(0, max_causal+1):
        cred_fname = f'{signal}/finemap_output.cred{n}'
        if os.path.exists(cred_fname):
            found_cred_file = True
            with open(cred_fname) as cred_file:
                first_line = next(cred_file)
            total_causal_count[n-1] += float(first_line.split()[-1])
            if n == max_causal:
                max_causal_loci.add(signal.split('/')[-1])
    if not found_cred_file:
        print(f'Failed to find cred_file for {signal.split("/")[-1]}')
        assert False

avg_causal_count = total_causal_count/n_signals_with_strs
assert np.abs(np.sum(avg_causal_count) - 1) < 0.01
fig, ax = plt.subplots()
ax.bar(np.arange(max_causal) + 1, avg_causal_count)
ax.set_xlabel('Number of causal varaints')
ax.set_ylabel('Average probability per locus')
ax.set_xticks(np.arange(max_causal) + 1)
plt.savefig(f'{ukb}/finemapping/finemap_results/{phenotype}/summary/avg_causal_count.png')

'''
for clump in clumps:
    per_var_output_df = pd.read_csv(
        f'{clump_dir(clump)}/finemap_output.snp',
        sep='\t',
        usecols=(1, 3, 10), #rsid, position, prob
        header=0
    )
'''

with open(f'{ukb}/finemapping/finemap_results/{phenotype}/summary/summary.txt', 'w') as summary:
    summary.write(f'Loci with some probability of the maximal number of causal variants ({max_causal}):\n')
    for locus in sorted(max_causal_loci):
        summary.write(f'{locus}\n')
    summary.write('\n')
