#!/usr/bin/env python3

import argparse

import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('n_vars_to_choose', type=int)
parser.add_argument('seed', type=int)
parser.add_argument('n_samples', type=int)
parser.add_argument('snp_mac_file')
# bins are 'bin_0001<001', 'bin_01>001', 'bin_1>01', 'bin_5>1'
parser.add_argument('bin_weights', nargs=4, type=float)
parser.add_argument('bin_effects', nargs=4)
#parser.add_argument('snp_gwas')
args = parser.parse_args()

rng = np.random.default_rng(seed=args.seed)

snps = pl.read_csv(
    args.snp_mac_file, sep='\t'
).with_column(
    pl.when(pl.col('mac') <= args.n_samples).then(pl.col('mac')).otherwise(args.n_samples - pl.col('mac'))
).select([
    (pl.col('mac')/args.n_samples/2).alias('maf'),
    (pl.col('chrom').cast(int).cast(str) + '_' + 'SNP' + '_' + pl.col('pos').cast(str) + '_' + pl.col('ref') + '_' + pl.col('alt')).alias('varname')
]).filter(
    pl.col('maf') >= 0.0001
).with_column(
    pl.when(pl.col('maf') < 0.001).
       then(args.bin_weights[0]).
       when((0.01 > pl.col('maf')) & (pl.col('maf') >= 0.001)).
       then(args.bin_weights[1]).
       when((0.1 > pl.col('maf')) & (pl.col('maf') >= 0.01)).
       then(args.bin_weights[2]).
       when((0.5 >= pl.col('maf')) & (pl.col('maf') >= 0.1)).
       then(args.bin_weights[3]).alias('weight')
)

assert snps.filter(pl.col('weight').is_null()).shape[0] == 0, snps.filter(pl.col('weight').is_null())

weights = snps['weight'].to_numpy().flatten()

idxs = rng.choice(snps.shape[0], size=args.n_vars_to_choose, replace=False, p=weights/np.sum(weights))

#gwas = pl.read_csv(args.snp_gwas, sep='\t', null_values='NA')

causal_vars = []
causal_betas = []

if args.n_vars_to_choose == 1:
    idxs = [idxs]

for idx in idxs:
    chrom, _, pos, ref, alt = snps[int(idx), 'varname'].split('_')
    causal_vars.append(f'SNP_{pos}_{ref}_{alt}')
    maf = snps[int(idx), 'maf']
    assert maf <= .5, (chrom, pos, ref, alt, maf)
    if maf >= 0.1:
        maf_idx = 3
    elif maf >= 0.01:
        maf_idx = 2
    elif maf >= 0.001:
        maf_idx = 1
    else:
        maf_idx = 0

    with open(args.bin_effects[maf_idx]) as bin_effects:
        effects = np.array([float(line) for line in bin_effects if line.strip()])
        causal_betas.append(rng.choice(effects))
#    snp_match = gwas.filter((pl.col('#CHROM') == int(chrom)) & (pl.col('POS') == int(pos)) & (pl.col('REF') == ref) & (pl.col('ALT') == alt))
#    assert snp_match.shape[0] == 1, snp_match
#    causal_betas.append(snp_match[0, 'BETA'])

print('\t'.join(causal_vars))
print('\t'.join(str(beta) for beta in causal_betas))
