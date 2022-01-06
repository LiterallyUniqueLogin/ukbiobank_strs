#!/usr/bin/env python3

import argparse

import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('chrom')
parser.add_argument('susie_dir')
parser.add_argument('mfi_file')
args = parser.parse_args()

colnames = pl.scan_csv(
    f'{args.susie_dir}/colnames.txt',
    sep='\t',
    has_header=False,
    with_column_names=lambda _: ['var_name'],
).filter(
    pl.col('var_name').str.contains('^SNP')
).with_column(
    pl.col('var_name').str.extract('SNP_(.*)')
).with_columns([
    pl.col('var_name').str.extract(r'^(\d+)_').cast(int).alias('pos'),
    pl.col('var_name').str.extract(r'^\d+_(\w+)_').alias('ref'),
    pl.col('var_name').str.extract(r'^\d+_\w+_(\w+)').alias('alt'),
])

mfi = pl.scan_csv(
    args.mfi_file,
    sep='\t',
    has_header=False,
    null_values='NA',
    with_column_names=lambda _: ['var_name', 'rsid', 'pos', 'ref', 'alt', 'mfi', 'alt2', 'info']
).filter(
    pl.col('rsid').str.contains('^rs')
)

use_vars = colnames.join(
    mfi,
    how='inner',
    on=['pos', 'ref', 'alt'],
).with_column(
    (1 - pl.col('mfi')).alias('opp_mfi')
).with_column(
    pl.when(
        pl.col('mfi') > 0.5
    ).then(
        pl.col('opp_mfi')
    ).otherwise(
        pl.col('mfi')
    ).alias('mfi')
).filter(
    (pl.col('mfi') > .05) & (pl.col('mfi') < .3)
).sort('pos').collect().to_dict()['rsid']#.sort('mfi', reverse=True).collect().head(150).to_dict()['rsid']

with open('out.txt', 'w') as out :
    for var in use_vars[::int(len(use_vars)/149)]:
        out.write(f'{var}\n')
    #[out.write(f'{var}\n') for var in use_vars]
    #[out.write(f'chr{args.chrom}:{var}\n') for var in use_vars]
