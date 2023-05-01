#!/usr/bin/env python3

import argparse

import bokeh.io
import bokeh.plotting
import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outdir')
parser.add_argument('locuscompare_tsv')
args = parser.parse_args()

df = pl.read_csv(
    args.locuscompare_tsv,
    sep='\t'
)

for col, step, bin_min, bin_max in zip(
    ['r', 'dosage_r', 'fraction_concordant_len', 'fraction_concordant_len_sum', 'numcalls', 'mean_absolute_difference'],
    [0.02, 0.02, 0.02, 0.02, 10000, 0.02],
    [None, None, None, None, None, 0],
    [1, 1, 1, 1, None, None],
):
#    if col == 'numcalls':
#        arr = df[col].to_numpy()
#    else:
#        arr = df.filter(pl.col('numcalls') > 150000)[col].to_numpy()
    arr = df[col].to_numpy()
    if bin_min is None:
        bin_min = np.min(np.floor(arr/step)*step)
    if bin_max is None:
        bin_max = np.max(np.ceil(arr/step)*step)
    bins = np.arange(bin_min, bin_max + step, step)
    fig = bokeh.plotting.figure(
        title=f'{col} histogram',
        y_axis_label='density',
        x_axis_label=col,
        width=1200,
        height=700,
        toolbar_location=None,
        output_backend='svg',
    )
    fig.grid.grid_line_color=None
    fig.quad(
        bottom=0,
        top=np.histogram(arr, bins=bins)[0]/len(arr),
        left=bins[:-1],
        right=bins[1:],
        line_width=0
    )
    bokeh.io.export_png(fig, filename=f'{args.outdir}/locuscompare_hist_{col}.png')
    bokeh.io.export_svg(fig, filename=f'{args.outdir}/locuscompare_hist_{col}.svg')

