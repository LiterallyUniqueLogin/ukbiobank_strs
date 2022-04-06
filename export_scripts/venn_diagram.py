#!/usr/bin/env python3

import argparse

import matplotlib.pyplot as plt
import matplotlib_venn
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outprefix')
parser.add_argument('results_tables', nargs='+')
args = parser.parse_args()

other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']

df = pl.concat([
    pl.scan_csv(
        table,
        sep='\t',
        dtypes={
                **{f'{ethnicity}_p_val': float for ethnicity in other_ethnicities},
                **{f'{ethnicity}_coeff': float for ethnicity in other_ethnicities},
                **{f'{ethnicity}_se': float for ethnicity in other_ethnicities}
        }
    ) for table in args.results_tables
]).filter(pl.col('p_val') <= 1e-10).with_columns([
    ((pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 0)).alias('susie_result'),
    (pl.col('finemap_pip') >= 0.8).alias('finemap_result')
]).filter(
    pl.col('susie_result') | pl.col('finemap_result')
).collect()

for var, condition in (('STR', pl.col('is_STR')), ('SNP', ~pl.col('is_STR'))):
    temp_df = df.filter(condition)
    s_total = temp_df.filter('susie_result').shape[0]
    f_total = temp_df.filter('finemap_result').shape[0]
    shared = temp_df.filter(pl.col('finemap_result') & pl.col('susie_result')).shape[0]

    plt.figure()
    matplotlib_venn.venn2(
        subsets = (s_total - shared, f_total - shared, shared),
        set_labels = ('SuSiE', 'FINEMAP')
    )

    plt.savefig(f'{args.outprefix}_{var}.png')
    plt.savefig(f'{args.outprefix}_{var}.svg')
