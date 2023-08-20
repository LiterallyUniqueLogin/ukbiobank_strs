#!/usr/bin/env python3

import argparse

import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('all_regions_finemapping_df')
parser.add_argument('n_samples', type=int)
parser.add_argument('mac_files', nargs='+')
args = parser.parse_args()

snps = pl.scan_csv(
    args.all_regions_finemapping_df, sep='\t'
).filter(
    ~pl.col('is_STR')
).select([
    (pl.col('chrom') + '_' + pl.col('varname')).alias('varname'),
    ((pl.col('finemap_pip') >= 0.5) | (
        (pl.col('susie_alpha') >= 0.5) &
        (pl.col('susie_cs') >= 0)
    )).alias('causal'),
    'coeff',
]).collect()

mafs = pl.concat([
    pl.read_csv(
        mac_file, sep='\t'
    ) for mac_file in args.mac_files
]).with_column(
    pl.when(pl.col('mac') <= args.n_samples).then(pl.col('mac')).otherwise(args.n_samples - pl.col('mac'))
).select([
    (pl.col('mac')/args.n_samples/2).alias('maf'),
    'mac',
    # fix chrom being '01' to just '1'
    (pl.col('chrom').cast(int).cast(str) + '_' + 'SNP' + '_' + pl.col('pos').cast(str) + '_' + pl.col('ref') + '_' + pl.col('alt')).alias('varname')
])

assert mafs.filter(pl.col('maf') > 0.5).shape[0] == 0, mafs.filter(pl.col('maf') > 0.5).shape[0]

# we're merging ignoring chromosome, so in theory
snps = snps.join(
    mafs,
    how='left',
    on='varname'
)

assert snps.filter(pl.col('maf').is_null()).shape[0] == 0, snps.filter(pl.col('maf').is_null())
assert snps.groupby('varname').agg(pl.count()).filter(pl.col('count') > 1).shape[0] == 0, snps.groupby('varname').agg(pl.count()).filter(pl.col('count') > 1)

snps = snps.filter(pl.col('maf') >= 0.0001)

snps = snps.with_columns([
    ((0.001 > pl.col('maf')) & (pl.col('maf') >= 0.0001)).alias('bin_001>0001'),
    ((0.01 > pl.col('maf')) & (pl.col('maf') >= 0.001)).alias('bin_01>001'),
    ((0.1 > pl.col('maf')) & (pl.col('maf') >= 0.01)).alias('bin_1>01'),
    ((0.5 >= pl.col('maf')) & (pl.col('maf') >= 0.1)).alias('bin_5>1'),
])

print('bin', 'num.', '% of total', 'num. causal', '% of total causal', 'per variant weight', sep='\t')
for bin_ in 'bin_001>0001', 'bin_01>001', 'bin_1>01', 'bin_5>1':
    print(
        bin_,
        snps.filter(bin_).shape[0],
        snps.filter(bin_).shape[0]/snps.shape[0],
        snps.filter(pl.col(bin_) & pl.col('causal')).shape[0],
        snps.filter(pl.col(bin_) & pl.col('causal')).shape[0]/snps.filter('causal').shape[0],
        snps.filter(pl.col(bin_) & pl.col('causal')).shape[0]/snps.filter('causal').shape[0]/snps.filter(bin_).shape[0],
        sep='\t'
    )
    with open(f'{bin_}_effects.txt'.replace('<', 'lthan').replace('>', 'gthan'), 'w') as effects:
        for coeff in snps.filter(pl.col(bin_) & pl.col('causal'))['coeff']:
            effects.write(f'{coeff}\n')
