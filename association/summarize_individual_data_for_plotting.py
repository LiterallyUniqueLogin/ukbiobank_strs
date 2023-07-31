#!/usr/bin/env python3

import argparse

import polars as pl
import statsmodels.stats.weightstats

parser = argparse.ArgumentParser()
parser.add_argument('out')
parser.add_argument('tab_file')
parser.add_argument('length_sum_column')
parser.add_argument('trait_column')
parser.add_argument('--sep', default='\t')
parser.add_argument('--subset', nargs=2)
args = parser.parse_args()

out = pl.read_csv(
    args.tab_file,
    sep=args.sep
)

if args.subset:
    out = out.filter(pl.col(args.subset[0]) == args.subset[1])

out = out.rename({
    args.length_sum_column: 'length_sum'
}).with_column(
    pl.col('length_sum').round(4)
).groupby(
    'length_sum'
).agg([
    pl.col(args.trait_column).mean().alias('mean'),
    pl.col(args.trait_column).apply(lambda vals: str(statsmodels.stats.weightstats.DescrStatsW(vals).tconfint_mean())).alias('conf_int'),
    pl.col(args.trait_column).count().alias('count')
]).sort('length_sum')

out.write_csv(args.out, sep='\t')

