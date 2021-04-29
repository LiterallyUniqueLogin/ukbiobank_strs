#!/usr/bin/env python3

import argparse
import os

import numpy as np

import load_and_filter_genotypes as lfg

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('phenotype')
parser.add_argument('chr')
parser.add_argument('--imputation-run-name')
parser.add_argument('--STRs', nargs='+')
parser.add_argument('--imputed-SNPs', nargs='+')

args = parser.parse_args()

assert len(args.STRs) + len(args.imputed_SNPs) > 0

if len(args.STRs) > 0:
    assert args.imputation_run_name

bgen_samples = []
with open(f'{ukb}/microarray/ukb46122_hap_chr1_v2_s487314.sample') as samplefile:
    for num, line in enumerate(samplefile):
        if num <= 1:
            # skip first two lines
            continue
        bgen_samples.append(line.split()[0])
assert len(bgen_samples) == 487409
samples_array = np.array(bgen_samples, dtype=float)

variant_names = []
gts = []

for STR in args.STRs:
    itr = lfg.load_strs(args.imputation_run_name, f'{args.chr}:{STR}-{STR}', slice(None))
    next(itr) # ignore names of locus_details
    try:
        dosage_dict, unique_alleles, chrom, pos, locus_filtered, locus_details = next(itr)
    except StopIteration as si:
        raise ValueError(f"Did not find a STR at position {STR}") from si

    assert pos == int(STR)
    if locus_filtered:
        raise ValueError(f"STR at position {STR} was filtered for reason {locus_filtered}")

    try:
        next(itr)
    except StopIteration:
        pass
    else:
        raise ValueError(f"Expected only one STR at position {STR}")

    variant_names.append(f"STR_{STR}")
    gts.append(np.sum([_len*np.sum(dosages, axis=1) for
                       _len, dosages in dosage_dict.items()], axis=0))

for iSNP in args.imputed_SNPs:
    itr = lfg.load_imputed_snps(f'{args.chr}:{iSNP}-{iSNP}', slice(None))
    next(itr) # ignore names of locus_details
    try:
        dosages, unique_alleles, chrom, pos, locus_filtered, locus_details = next(itr)
    except StopIteration as si:
        raise ValueError(f"Did not find an imputed SNP at position {iSNP}") from si

    assert pos == int(iSNP)
    if locus_filtered:
        raise ValueError(f"imputed SNP at position {iSNP} was filtered for reason {locus_filtered}")

    try:
        next(itr)
    except StopIteration:
        pass
    else:
        raise ValueError(f"Expected only one imputed SNP at position {iSNP}")

    variant_names.append(f"imputed_SNP_{iSNP}")
    gts.append(dosages[:, 1] + 2*dosages[:, 2])

STR_string = '_'.join(args.STRs)
if len(STR_string) > 0:
    STR_string = '_' + STR_string

iSNP_string = '_'.join(args.imputed_SNPs)
if len(iSNP_string) > 0:
    iSNP_string = '_' + iSNP_string

fname = f'chr{args.chr}_STR{STR_string}__ISNP{iSNP_string}__ASNP'

gts.insert(0, samples_array)
dirname = f'{ukb}/association/results/height/conditional_inputs'
np.save(
    f'{dirname}/{fname}.npy',
    np.stack(gts, axis=1)
)
with open(f'{dirname}/{fname}_varnames.txt', 'w') as f:
    f.write('sample_ID ' + ' '.join(variant_names) + '\n')

with open(f'{dirname}/{fname}_README.txt', 'w') as f:
    f.write('Computed single dosages per sample for STRs in number of repeats, single '
            'dosages per sample for SNPs as number of copies of the alternate allele.\n')
