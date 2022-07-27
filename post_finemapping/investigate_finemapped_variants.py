#!/usr/bin/env python3

import argparse
import os

import polars as pl

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('name')
parser.add_argument('--phenotypes', nargs='+')
parser.add_argument('--regions', nargs='+')
args = parser.parse_args()

assert len(args.phenotypes) == len(args.regions)
dfs =  [pl.read_csv(
    f'{ukb}/post_finemapping/intermediate_results/finemapping_all_concordance_{pheno}.tab',
    sep='\t'
).filter(
    (pl.col('p_val') <= 5e-8) &
    (pl.col('region') == region)
).filter(
    (pl.col('finemap_pip') > 0.01) | (
        (pl.col('susie_alpha') > 0.01) & (pl.col('susie_cs') > 0)
    )
).select([
    'varname',
    'p_val',
    'susie_cs',
    'susie_alpha',
    'finemap_pip'
]) for pheno, region in zip(args.phenotypes, args.regions)]

joined_df = dfs[0]
joined_df = joined_df.rename({colname: f'{colname}_{args.phenotypes[0]}' for colname in ['p_val', 'susie_cs', 'susie_alpha', 'finemap_pip']})
for df_count, df in enumerate(dfs[1:]):
    df = df.rename({colname: f'{colname}_{args.phenotypes[df_count + 1]}' for colname in ['p_val', 'susie_cs', 'susie_alpha', 'finemap_pip']})
    joined_df = joined_df.join(
        df,
        how='outer',
        on='varname'
    )

joined_df.write_csv(f'{args.name}_finemapped_variants.tab', sep='\t')
