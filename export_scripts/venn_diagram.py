#!/usr/bin/env python3

import argparse

import matplotlib.pyplot as plt
import matplotlib_venn
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outprefix')
parser.add_argument('results_table')
args = parser.parse_args()

df = pl.scan_csv(
    args.results_table,
    sep='\t'
).filter(pl.col('p_val') <= 1e-10).with_columns([
    ((pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 0)).alias('susie_result'),
    (pl.col('finemap_pip') >= 0.8).alias('finemap_result')
]).filter(
    pl.col('susie_result') | pl.col('finemap_result')
).collect()

s_total = df.filter('susie_result').shape[0]
f_total = df.filter('finemap_result').shape[0]
shared = df.filter(pl.col('finemap_result') & pl.col('susie_result')).shape[0]

matplotlib_venn.venn2(
    subsets = (s_total - shared, f_total - shared, shared),
    set_labels = ('SuSiE', 'FINEMAP')
)

plt.savefig(f'{args.outprefix}.png')
plt.savefig(f'{args.outprefix}.svg')
