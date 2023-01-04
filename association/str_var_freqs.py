#!/usr/bin/env python3

import argparse

import cyvcf2
import numpy as np

import sample_utils

parser = argparse.ArgumentParser()
parser.add_argument('str_vcf')
parser.add_argument('chrom', type=int)
parser.add_argument('pos', type=int)
parser.add_argument('all_samples_fname')
# in the order white brits, black, south asian, chinese
parser.add_argument('ethnic_samples_fnames', nargs=4)

args = parser.parse_args()

var = next(cyvcf2.VCF(args.str_vcf)(f'{args.chrom}:{args.pos}'))

ethnicities = ('white_brits', 'black', 'south_asian', 'chinese')
print([var.REF] + var.ALT)
for ethnicity, ethnic_samples_fname in zip(ethnicities, args.ethnic_samples_fnames):
    samp_idx = sample_utils.get_samples_idx(args.all_samples_fname, ethnic_samples_fname)
    print(f'n_samples {ethnicity} = {np.sum(samp_idx)}')

print('ethnicity\tallele\tfreq')
for ethnicity, ethnic_samples_fname in zip(ethnicities, args.ethnic_samples_fnames):
    samp_idx = sample_utils.get_samples_idx(args.all_samples_fname, ethnic_samples_fname)

    freq = np.sum((
        np.maximum(0, 1 - np.sum(var.format('AP1'), axis=1)) +
        np.maximum(0, 1 - np.sum(var.format('AP2'), axis=1))
    )[samp_idx])/(2*np.sum(samp_idx))
    print(f'{ethnicity}\t0\t{freq:.3}')

    for idx, seq in enumerate(var.ALT):
    	freq = np.sum((var.format('AP1') + var.format('AP2'))[samp_idx, idx])/(2*np.sum(samp_idx))
    	print(f'{ethnicity}\t{idx+1}\t{freq:.3}')

