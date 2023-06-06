#!/usr/bin/env python3

import argparse

import bgen_reader
import numpy as np

import load_and_filter_genotypes as lfg
import sample_utils

parser = argparse.ArgumentParser()
parser.add_argument('out')
parser.add_argument('str_vcf')
parser.add_argument('snp_bgen')
parser.add_argument('all_samples_list')
parser.add_argument('samples_list')
parser.add_argument('chrom')
parser.add_argument('seed', type=int)
parser.add_argument('--causal-vars', nargs='*') # either POS which is an STR or POS_REF_ALT which is a SNP
parser.add_argument('--causal-betas', nargs='*', type=float)
args = parser.parse_args()

assert len(args.causal_vars) == len(args.causal_betas)

sample_idx = sample_utils.get_samples_idx(args.all_samples_list, args.samples_list)
samples = sample_utils.get_ordered_samples(args.all_samples_list, args.samples_list)
assert len(samples) == np.sum(sample_idx)

explained_variance = 0
phenotypes = np.zeros(np.sum(sample_idx))

bgen = bgen_reader.open_bgen(args.snp_bgen, verbose=False)

for var, beta in zip(args.causal_vars, args.causal_betas):
    print(var, beta)
    split = var.split('_')
    pos = int(split[1])
    if split[0] == 'STR':
        itr = lfg.load_strs(
            args.str_vcf,
            f'{args.chrom}:{pos}-{pos}',
            sample_idx,
            details=False
        )
        probs, len_alleles, _, __, ___, _____ = next(itr)
        caught = False
        try:
            next(itr)
        except StopIteration:
            caught = True
        assert caught

        dosages = sum(len_ * np.sum(len_prob, axis=1) for len_, len_prob in probs.items())
        dosages = dosages - np.mean(dosages)

        phenotypes += beta*dosages
        explained_variance += np.std(dosages)**2 * beta**2
    else:
        assert split[0] == 'SNP'

        ref, alt = split[2:]
        where = np.where((bgen.positions == pos) & (bgen.allele_ids == f'{ref},{alt}'))[0]
        assert len(where) == 1
        probs = bgen.read(where[0]).squeeze()[sample_idx, :]

        dosages = probs[:, 1] + 2*probs[:, 2]
        dosages = dosages - np.mean(dosages)

        phenotypes += beta*dosages
        explained_variance += np.std(dosages)**2 * beta**2

assert explained_variance < 0.1

rng = np.random.default_rng(args.seed)
phenotypes += rng.normal(0, np.sqrt(1 - explained_variance), phenotypes.shape)
phenotypes = (phenotypes - np.mean(phenotypes))/np.std(phenotypes)

samp_phen = np.hstack((samples.reshape(-1, 1), phenotypes.reshape(-1, 1)))
np.save(args.out, samp_phen)
