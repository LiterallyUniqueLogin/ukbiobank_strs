#!/usr/bin/env python3

import argparse

import bokeh.plotting
import numpy as np
import polars as pl
import scipy.stats

parser = argparse.ArgumentParser()
parser.add_argument('outdir')
parser.add_argument('chrom_files', nargs = '+',
                    help='4 cols: pos, chance of length confusionn, avg abs length confusion, normalized avg abs lenght confusion')
args = parser.parse_args()
outdir = args.outdir
chrom_fnames = args.chrom_files

loci = pl.concat([
    pl.scan_csv(
        chrom_fname,
        sep='\t'
    ) for chrom_fname in chrom_fnames
]).drop('pos').collect()

for col in loci.columns:
    print(f'Plotting column {col} ...', flush=True)
    max_val = loci.select(pl.col(col).max()).to_numpy()
    min_val = loci.select(pl.col(col).min()).to_numpy()
    n_steps = 1000
    step_size = (max_val - min_val)/n_steps
    xs = np.arange(min_val, max_val + step_size, step_size)
    ys = scipy.stats.gaussian_kde(loci[col].to_numpy())(xs)

    if col.startswith('chance'):
        unit = '%'
    elif col.startswith('avg'):
        unit = 'repeat units'
    else:
        assert col.startswith('normalized')
        unit = 'length standard deviations'
    x_axis_label = f'{col} ({unit})'
    
    figure = bokeh.plotting.figure(
        title=col.capitalize() + ' distribution',
        x_axis_label=x_axis_label,
        y_axis_label='pdf',
        width=1200,
        height=1200,
        tools=''
    )
    figure.background_fill_color = None
    figure.border_fill_color = None
    figure.grid.grid_line_color = None
    figure.title.text_font_size='30px'
    figure.axis.axis_label_text_font_size='26px'
    figure.axis.major_label_text_font_size='20px'

    figure.varea(
        x=xs,
        y1=np.zeros(xs.shape),
        y2=ys
    )

    bokeh.io.export_png(figure, filename=f'{outdir}/{col}_distribution.png')
    bokeh.io.export_svg(figure, filename=f'{outdir}/{col}_distribution.svg')

print('Done')
