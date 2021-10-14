#!/usr/bin/env python3

import glob
import os

import numpy as np

ukb = os.environ['UKB']

def to_float(f):
    if f == 'nan':
        return np.nan
    else:
        return -np.log10(max(float(f), 1e-300))

snps = []
strs = []
fnames = glob.glob(f'{ukb}/signals/peaks/*_250000_5e-8.tab')
for fname in fnames:
    with open(fname) as peaks:
        header = next(peaks).strip().split('\t')
        var_type_idx = header.index('variant_type')
        p_val_idx = header.index('p_value')
        p_val_other_idx = header.index('p_value_other_variant_type')
        for line in peaks:
            var_type, p_val, other_p_val = np.array(line.strip().split())[
                [var_type_idx, p_val_idx, p_val_other_idx]
            ]
            p_val = to_float(p_val)
            other_p_val = to_float(other_p_val)

            if np.isnan(p_val):
                p_val = 0
            if np.isnan(other_p_val):
                other_p_val = 0

            if var_type == 'SNP':
                snps.append(p_val)
                strs.append(other_p_val)
            else:
                assert var_type == 'STR'
                snps.append(other_p_val)
                strs.append(p_val)
snps = np.array(snps)
strs = np.array(strs)

print(f'n peaks: {len(snps)}')
print(f'n peaks per pheno: {len(snps)/len(fnames)}')
print(f'% peaks only STR: {np.sum(snps == 0)/len(snps)}')

print('Subsetting to peaks that have a SNP')
strs = strs[snps != 0]
snps = snps[snps != 0]

for cutoff in .8, .9, 1:
    print(f'fraction STR peaks >= {cutoff} fraction of SNP peak {np.sum(snps*cutoff <= strs)/len(snps)}')
    for sig_thresh in 20, 30, 50:
        sig_idx = (snps >= sig_thresh) | (strs >= sig_thresh)
        print(f'fraction STR peaks >= {cutoff} fraction of SNP peak among peaks <= 1e-{sig_thresh} {np.sum((snps*cutoff <= strs)[sig_idx])/np.sum(sig_idx)}')


non_capped_idx = (snps < 300) & (strs < 300)

snps = snps[non_capped_idx]
strs = strs[non_capped_idx]

for cutoff in .8, .9, 1:
    print(f'fraction STR peaks >= {cutoff} fraction of SNP peak {np.sum(snps*cutoff <= strs)/len(snps)}')
    for sig_thresh in 20, 30, 50:
        sig_idx = (snps >= sig_thresh) | (strs >= sig_thresh)
        print(f'fraction STR peaks >= {cutoff} fraction of SNP peak among peaks <= 1e-{sig_thresh} {np.sum((snps*cutoff <= strs)[sig_idx])/np.sum(sig_idx)}')

