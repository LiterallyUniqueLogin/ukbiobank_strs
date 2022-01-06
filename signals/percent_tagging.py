#!/usr/bin/env python3

import os

import numpy as np

import phenotypes

ukb = os.environ['UKB']

def to_float(f):
    if f == 'nan':
        return np.nan
    else:
        return -np.log10(max(float(f), 1e-300))

snps = []
strs = []
fnames = [f'{ukb}/signals/peaks/{phenotype}_250000_5e-8.tab' for phenotype in phenotypes.phenotypes_in_use]
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
                p_val = 1
            if np.isnan(other_p_val):
                other_p_val = 1

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
print(f'% peaks only STR: {np.sum(snps == 1)/len(snps)}')
print(f'% peaks only SNP: {np.sum(strs == 1)/len(strs)}')

print('Strongest STR peak with no SNP:', np.max(strs[snps == 1]))

