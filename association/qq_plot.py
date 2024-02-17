#!/usr/bin/env python3

import argparse

import bokeh.io
import bokeh.plotting
import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('results_tab')
parser.add_argument('p_val_col')
parser.add_argument('type')
parser.add_argument('pheno')
parser.add_argument('out')
parser.add_argument('--null-values')
parser.add_argument('--y-max', type=float)
args = parser.parse_args()

p_vals = pl.scan_csv(
    args.results_tab,
    separator='\t',
    null_values = args.null_values
).select(
    pl.col(args.p_val_col).alias('p_val')
).filter(
    ~pl.col('p_val').is_null()
).select(
    (-pl.when(
        pl.col('p_val') < 1e-300
    ).then(
        1e-300
    ).otherwise(
        pl.col('p_val')
    ).log10()).alias('p_val')
).sort('p_val').collect().to_numpy().flatten()

if args.y_max:
    p_vals = np.minimum(p_vals, args.y_max)

null_dist_vals = -np.log10(
    np.arange(0, 1, 1/len(p_vals)) + 1/(len(p_vals)*2)
)[::-1]

fig = bokeh.plotting.figure(
    x_range=bokeh.models.Range1d(0, max(null_dist_vals)),
    y_range=bokeh.models.Range1d(0, max(p_vals) if args.y_max is None else args.y_max),
    title=f'{args.type} {args.pheno} QQ plot',
    tools='',
    x_axis_label = 'uniform distribution p-values',
    y_axis_label = 'measured p-values',
    width=1200,
    height=1200,
)
line_max = min(max(null_dist_vals), max(p_vals))
fig.line(
    x=[0, line_max],
    y=[0, line_max]
)
fig.circle(
    x=null_dist_vals[::10000],
    y=p_vals[::10000]
)
fig.circle(
    x=null_dist_vals[-500_000:],
    y=p_vals[-500_000:]
)
bokeh.io.export_png(fig, filename=f'{args.out}.png')
