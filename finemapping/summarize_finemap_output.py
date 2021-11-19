#!/usr/bin/env python3

import argparse
import os
import os.path

import matplotlib.pyplot as plt
import numpy as np

import graphing_utils

parser = argparse.ArgumentParser()
parser.add_argument('phenotype')
parser.add_argument('out_dir')
parser.add_argument('n_signals_no_strs')
parser.add_argument('signal_dirs', nargs='+')
args = parser.parse_args()
phenotype = args.phenotype
n_signals_no_strs = int(args.n_signals_no_strs)
signal_dirs = args.signal_dirs
out_dir = args.out_dir
n_signals_with_strs = len(signal_dirs)
n_signals = n_signals_with_strs + n_signals_no_strs

# ---- what do I want to calculate
# ranking of most significant STRs across all loci
# percent and raw count loci with an STR in the most probable configuration, and list of such
# percent and raw count STRs that reach a certain threshold (ROC curve)
# percent and raw count SNPs that reach a certain threshold (ROC curve)
# loci ranked by highest probability - call out loci with low prob
# loci with number of causal SNPs - loci with too many causal SNPs
# reproducability?

# number of causal variants per locus

max_causal = 20

total_causal_count = np.zeros((max_causal+1))

max_causal_loci = set()

for signal_dir in signal_dirs:
    found_cred_file = False
    for n in range(0, max_causal+1):
        cred_fname = f'{signal_dir}/finemap_output.cred{n}'
        if os.path.exists(cred_fname):
            found_cred_file = True
            with open(cred_fname) as cred_file:
                first_line = next(cred_file)
            total_causal_count[n] += float(first_line.split()[-1])
            if n == max_causal:
                max_causal_loci.add(signal_dir.split('/')[-1])
    if not found_cred_file:
        print(f'Failed to find cred_file for {signal_dir.split("/")[-1]}')
        assert False

avg_causal_count = total_causal_count/np.sum(total_causal_count)
fig, ax = plt.subplots()
ax.bar(np.arange(max_causal + 1), avg_causal_count)
ax.set_xlabel('Number of causal varaints')
ax.set_ylabel('Average probability per locus')
ax.set_xticks(np.arange(max_causal + 1))
plt.savefig(f'{out_dir}/avg_causal_count.png')

str_pcausals = []
pcausals = []
highest_single_str_pcausal = []
highest_str_rank = []
highest_str_id = []

all_str_pcausal = []
all_str_rank = []
all_str_id = []
all_str_signal = []
for signal_dir in signal_dirs:
    signal_name = signal_dir.split('/')[-1]
    pcausals.append(0)
    str_pcausals.append(0)
    any_strs = False
    with open(f'{signal_dir}/finemap_output.snp') as per_var_output:
        next(per_var_output)
        for rank, line in enumerate(per_var_output):
            rank += 1
            splits = line.split()
            id_ = splits[1]
            prob = float(splits[10])
            pcausals[-1]  += prob
            if id_[:4] == 'SNP_':
                continue
            elif id_[:4] != 'STR_':
                raise ValueError(f'Unidentified id {id_}')
            str_pcausals[-1] += prob
            if not any_strs:
                any_strs = True
                highest_single_str_pcausal.append(prob)
                highest_str_rank.append(rank)
                highest_str_id.append(id_)
            all_str_pcausal.append(prob)
            all_str_id.append(id_)
            all_str_rank.append(rank)
            all_str_signal.append(signal_name)
        if not any_strs:
            highest_single_str_pcausal.append(0)
            highest_str_rank.append(np.inf)
            highest_str_id.append('NA')

total_pcausal = np.sum(pcausals)
total_str_pcausal = np.sum(str_pcausals)
str_contrib_fraction = total_str_pcausal/total_pcausal
str_contrib_fractions = [
    str_contrib/total_pcausal for (str_contrib, total_pcausal) in zip(str_pcausals, pcausals)
]

graphing_utils.PlotCDF(
    np.array(str_contrib_fractions),
    f'{out_dir}/str_contrib_fractions.png',
    'signals',
    'fractional STR contribution to causality',
    data_len = n_signals
)
graphing_utils.PlotCDF(
    np.array(highest_single_str_pcausal),
    f'{out_dir}/single_str_contrib_fractions.png',
    'signals',
    "largest single STR causal probability",
    data_len = n_signals
)
graphing_utils.PlotCDF(
    np.array(highest_str_rank),
    f'{out_dir}/str_rank.png',
    'signals',
    "rank of most probable STR",
    reverse = True
)

with open(f'{out_dir}/best_STR_contribs.tab', 'w') as best_contribs:
    best_contribs.write('The best STR from each signal, sorted by their posterior probability of causality\n')
    best_contribs.write('signal\tSTR\tpcausal\trank\n')
    best_contrib_order = np.argsort(highest_single_str_pcausal, kind='stable')[::-1]
    for idx in best_contrib_order:
        pcausal = highest_single_str_pcausal[idx]
        rank = highest_str_rank[idx]
        id_ = highest_str_id[idx]
        signal = signal_dirs[idx].split('/')[-1]
        best_contribs.write(f'{signal}\t{id_}\t{pcausal}\t{rank}\n')

with open(f'{out_dir}/best_STR_ranks.tab', 'w') as best_ranks:
    best_ranks.write('The best STR from each signal, sorted by their ranks\n')
    best_ranks.write('(the rank of an STR the rank of its causal probability compared to other variants at the signal of that STR)\n')
    best_ranks.write('signal\tSTR\trank\tpcausal\n')
    best_contrib_order = np.argsort(highest_str_rank, kind='stable')
    for idx in best_contrib_order:
        pcausal = highest_single_str_pcausal[idx]
        rank = highest_str_rank[idx]
        id_ = highest_str_id[idx]
        signal = signal_dirs[idx].split('/')[-1]
        best_ranks.write(f'{signal}\t{id_}\t{rank}\t{pcausal}\n')

with open(f'{out_dir}/all_STR_contribs.tab', 'w') as all_contribs:
    all_contribs.write('All STRs from all signals sorted by posterior probability of causality\n')
    #all_contribs.write('Cutoff is pcausal >= 0.05\n')
    all_contribs.write('signal\tSTR\tpcausal\trank\n')
    all_str_order = np.argsort(all_str_pcausal, kind='stable')[::-1]
    for idx in all_str_order:
        pcausal = all_str_pcausal[idx]
        id_ = all_str_id[idx]
        rank = all_str_rank[idx]
        signal = all_str_signal[idx]
        all_contribs.write(f'{signal}\t{id_}\t{pcausal}\t{rank}\n')

with open(f'{out_dir}/summary.txt', 'w') as summary:
    summary.write(
        f'{n_signals_no_strs} signals without str variants at all above threshold 1e-4, {n_signals_with_strs} signals with such STR variants, '
        f'{n_signals_no_strs + n_signals_with_strs} total variants.\n\n'
    )
    if len(max_causal_loci) == 0:
        summary.write(f'No loci with any probability of the maximal number of causal variants ({max_causal})\n')
    else:
        summary.write(f'Loci with some probability of the maximal number of causal variants ({max_causal}):\n')
        for locus in sorted(max_causal_loci):
            summary.write(f'{locus}\n')
    summary.write('\n')

    summary.write(
        'Percent causality contributed by STRs (no filtering for low probability variants): '
        f'{str_contrib_fraction*100:.2f}%'
    )
