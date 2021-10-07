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

tag_fractions = []
fnames = glob.glob(f'{ukb}/signals/peaks/*_250000_5e-8.tab')
for fname in fnames:
    max_p = 0
    is_STR = False
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

            if p_val == other_p_val == 300:
                min_p = np.nan
            elif p_val > max_p:
                max_p = p_val
                is_STR = var_type == 'STR'

            if var_type == 'SNP':
                if np.isnan(other_p_val):
                    tag_fractions.append(0)
                else:
                    tag_fractions.append(other_p_val/p_val)
            else:
                assert var_type == 'STR'
                if np.isnan(other_p_val):
                    tag_fractions.append(np.inf)
                else:
                    tag_fractions.append(p_val/other_p_val)
    if is_STR:
        print(fname)


print(f'n peaks: {len(tag_fractions)}')
print(f'peaks per pheno: {len(tag_fractions)/len(fnames)}')

tag_fractions = np.array(tag_fractions)

print(f'fraction of peaks with only STRs')
print(np.sum(np.isinf(tag_fractions))/tag_fractions.shape[0])

for cutoff in .8, .9, 1:
    print(f'fraction bigger than {cutoff}')
    print(np.sum(tag_fractions >= cutoff)/tag_fractions.shape[0])
