#!/usr/bin/env python3

import argparse

import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('str_stats')
parser.add_argument('snp_stats')
parser.add_argument('snp_filter_files')
args = parser.parse_args()

snps_to_filter_df = pl.read_csv(
    args.snp_filter_files,
    sep='\t'
).select([
    'snp_pos', 'snp_ref', 'snp_alt'
])
snps_to_filter = list(zip(*(snps_to_filter_df[col] for col in snps_to_filter_df.columns)))

snps = pl.read_csv(
    args.snp_stats,
    sep='\t'
)
haplotype_count = snps[0, 'OBS_CT']
print('haplotype_count', haplotype_count)

snp_nonmajor_allele_freq = 0
snp_nonmajor_basepair_freq = 0
snp_count = 0
for i in range(snps.shape[0]):
    if (snps[i, 'POS'], snps[i, 'REF'], snps[i, 'ALT']) in snps_to_filter:
        continue

    snp_count += 1

    freq = min(snps[i, 'ALT_CTS']/haplotype_count, 1 - snps[i, 'ALT_CTS']/haplotype_count)
    snp_nonmajor_allele_freq += freq

    ref = snps[i, 'REF']
    alt = snps[i, 'ALT']

    shared_front = 0
    while shared_front < len(ref) and shared_front < len(alt) and ref[shared_front] == alt[shared_front]:
        shared_front += 1
    if shared_front > 0:
        ref = ref[shared_front:]
        alt = alt[shared_front:]

    shared_back = 0
    while -shared_back < len(ref) and -shared_back < len(alt) and ref[shared_back-1] == alt[shared_back-1]:
        shared_back -= 1
    if shared_back < 0:
        ref = ref[:shared_back]
        alt = alt[:shared_back]

    # this is a SNP
    if len(ref) == 1 and len(alt) == 1:
        snp_nonmajor_basepair_freq += freq
        continue

    # this is an indel
    assert len(ref) == 0 or len(alt) == 0, (snps[i, 'POS'], snps[i, 'REF'], snps[i, 'ALT'], ref, alt)

    snp_nonmajor_basepair_freq += freq*max(len(ref), len(alt))

print('snp_count', snp_count)
print('snp_mean_nonmajor_alleles_per_haplotype', snp_nonmajor_allele_freq)
print('snp_mean_nonmajor_basepairs_per_haplotype', snp_nonmajor_basepair_freq, flush=True)

str_nonmajor_allele_freq = 0
str_nonmajor_basepair_freq = 0
strs = pl.read_csv(args.str_stats, sep='\t').with_columns([
    pl.col('alleles').str.split(',').cast(pl.List(int)),
    pl.col('counts').str.split(',').cast(pl.List(float)),
])
assert haplotype_count == strs[0, 'n']*2, (haplotype_count, strs[0, 'n']*2)
for i in range(strs.shape[0]):
    alleles = strs[i, 'alleles'].to_numpy()
    counts = strs[i, 'counts'].to_numpy()
    max_idx = np.argmax(counts)
    str_nonmajor_allele_freq += 1 - counts[max_idx]/haplotype_count
    str_nonmajor_basepair_freq += np.sum(np.abs(alleles - alleles[max_idx])*counts)/haplotype_count

print('str_count', strs.shape[0])
print('str_mean_nonmajor_allele_per_haplotype', str_nonmajor_allele_freq)
print('str_mean_nonmajor_basepair_per_haplotype', str_nonmajor_basepair_freq, flush=True)

