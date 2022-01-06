#!/usr/bin/env python3

# Solely based on dosages
# ROC for each allele
# AUC for each allele
# confusion matrix: weighted by dosages

import argparse
import collections
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import load_and_filter_genotypes as lfg
import python_array_utils as utils

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('chrom')
parser.add_argument('pos')
parser.add_argument('samples_fname')
parser.add_argument('outfname')
args = parser.parse_args()

with open(args.samples_fname) as samples_file:
    samples = np.array([line.strip() for line in samples_file][1:], dtype=int).reshape(-1, 1)

imp_snp_samples_filepath = f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample'
with open(imp_snp_samples_filepath) as imp_snp_samples_file:
    imp_snp_samples = np.array([line.split()[0] for line in imp_snp_samples_file][2:], dtype=int).reshape(-1, 1)

samples_indicator = np.concatenate((samples, samples), axis=1)
samples_merge = utils.merge_arrays(imp_snp_samples, samples_indicator)
assert samples_merge.shape[1] == 2
sample_idx = ~np.isnan(samples_merge[:, 1])

strs = lfg.load_strs(
    'first_pass',
    f'{args.chrom}:{args.pos}-{args.pos}',
    sample_idx,
    details = False
)

single_chrom_dosages = next(strs)[0]
summed_dosages = collections.defaultdict(lambda: np.zeros(np.sum(sample_idx)))
for len1, dosages1 in single_chrom_dosages.items():
    #print(len1, dosages1)
    for len2, dosages2 in single_chrom_dosages.items():
        if len2 < len1:
            continue
        #print('------', len2, dosages2)
        summed_dosages[round(len1 + len2, 2)] += (
            np.multiply(dosages1[:, 0], dosages2[:, 1]) +
            np.multiply(dosages2[:, 0], dosages1[:, 1])
        )
#print(summed_dosages)

single_chrom_dosages = {
    k: np.hstack((v[:, 0], v[:, 1])) for (k, v) in single_chrom_dosages.items()
}

f, ax = plt.subplots(2, 2, figsize=(32, 32))
ax[0, 0].set_title('Confusion matix')
ax[0, 1].set_title('Normalized confusion matix')
for row, dosages in enumerate((single_chrom_dosages, summed_dosages)):
    # TODO make threshold more apparent if we want to display this
    dosages = {k: v for (k,v) in dosages.items() if np.sum(v) >= 100}
    alleles = sorted(dosages)
    n_alleles = len(alleles)

    confusion = np.full((n_alleles, n_alleles), np.nan)

    for i in range(n_alleles):
        dosage_i = dosages[alleles[i]]
        for j in range(n_alleles):
            dosage_j = dosages[alleles[j]]
            confusion[i, j] = np.dot(dosage_i, dosage_j)
            '''
            confusion[i, j] = (
                np.dot(dosage_i[:, 0], dosage_j[:, 0]) +
                np.dot(dosage_i[:, 1], dosage_j[:, 1])
            )
            '''

    assert not np.any(np.isnan(confusion))

    normalized_confusion = confusion.copy()
    for i in range(n_alleles):
        normalized_confusion[i, :] /= np.sum(normalized_confusion[i, :])

    confusion_df = pd.DataFrame(confusion, index=alleles, columns=alleles)
    confusion_df.axes[0].name = 'predicted allele'
    confusion_df.axes[1].name = 'true allele'
    normalized_confusion_df = pd.DataFrame(normalized_confusion, index=alleles, columns=alleles)
    normalized_confusion_df.axes[0].name = 'predicted allele'
    normalized_confusion_df.axes[1].name = 'true allele'

    sns.heatmap(confusion_df, annot=True, fmt='.0f', linewidths=.5, ax=ax[row, 0], square=True)

    sns.heatmap(normalized_confusion_df, vmin=0, vmax=1, annot=True, fmt=".2%", linewidths=.5, ax=ax[row, 1], square=True)
ax[0, 0].set_ylabel('Per Allele\npredicted allele')
ax[1, 0].set_ylabel('Summed Genotype\npredicted allele')

plt.savefig(args.outfname)
plt.savefig(args.outfname.replace('png', 'pdf'))
