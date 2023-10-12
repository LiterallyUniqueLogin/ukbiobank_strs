#!/usr/bin/env python3

import argparse

import bokeh.io
import bokeh.layouts
import bokeh.plotting
import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outdir')
parser.add_argument('locuscompare_tsvs', nargs=6)
parser.add_argument('subset_chr_pos')
args = parser.parse_args()

subset = pl.read_csv(
    args.subset_chr_pos,
    sep='\t'
).unique()

dfs = []
for locuscompare_tsv in args.locuscompare_tsvs:
    df = pl.read_csv(
        locuscompare_tsv,
        sep='\t'
    ).with_columns([
        (pl.col('r')**2).alias('r2'),
        (pl.col('dosage_r')**2).alias('dosage_r2')
    ])

    df = subset.join(
        df,
        how='inner',
        left_on=['chrom', 'start_pos (hg19)'],
        right_on=['chrom', 'start'],
    )

    assert df.shape[0] == subset.shape[0]
    assert df.shape[0] == 409
    dfs.append(df)

rows = []
nice_cols = ['Number of WGS calls', 'Fraction of calls with\nconcordant length sum', 'Mean absolute\nlength sum difference\n(repeat units)', 'dosage r^2'] # fraction_concordant_len, 'weighted r^2'
for df, ethnicity in zip(dfs, ['White British', 'Black', 'South Asian', 'Chinese', 'Irish', 'White Other']):
    row = []
    for idx, (col, nice_col, step, bin_min, bin_max) in enumerate(zip(
        ['numcalls', 'fraction_concordant_len_sum', 'mean_absolute_difference', 'dosage_r2'], #'r2',  fraction_concordant_len
        nice_cols,
        [None, 0.02, 0.10, 0.02], # 0.02, 0.02
        [None, None, 0, None], # None, None
        [None, 1, None, 1], # 1, 1
    )):
        arr = df[col].to_numpy()

        if col == 'dosage_r2':
            print(df.select([
                (pl.col(col) > 0.9).cast(int).mean().alias('0.9'),
                (pl.col(col) > 0.8).cast(int).mean().alias('0.8'),
                (pl.col(col) > 0.6).cast(int).mean().alias('0.6'),
            ]))

        if col == 'numcalls':
            if ethnicity != 'White British':
                step = 50
            else:
                step = 4000

        if bin_min is None:
            if col != 'numcalls':
                bin_min = min(np.min(np.floor(df_[col].to_numpy()/step)*step) for df_ in dfs)
            else:
                bin_min = np.min(np.floor(arr/step)*step)
        if bin_max is None:
            if col != 'numcalls':
                bin_max = max(np.max(np.ceil(df_[col].to_numpy()/step)*step) for df_ in dfs)
            else:
                bin_max = np.max(np.ceil(arr/step)*step)

        bins = np.arange(bin_min, bin_max + step, step)
        fig_kwargs = dict(
            width=1200,
            height=700,
            toolbar_location=None,
            output_backend='svg',
        )
        if idx == 0:
            fig_kwargs['y_axis_label'] = ethnicity
#        if len(rows) == 5:
#            fig_kwargs['x_axis_label'] = nice_col
        if col != 'numcalls':
            y_max = max(max(np.histogram(df_[col].to_numpy(), bins=bins)[0]/df_.shape[0]) for df_ in dfs)
            fig_kwargs['y_range'] = (0, 1.025*y_max)
        fig = bokeh.plotting.figure(
            **fig_kwargs
        )
        if idx == 0:
            fig.xaxis[0].formatter = bokeh.models.NumeralTickFormatter(format="0a")
        fig.grid.grid_line_color = None
        fig.axis.axis_label_text_font_size = '72px'
        fig.axis.major_label_text_font_size = '48px'
        fig.title.text_font_size = '72px'
        fig.quad(
            bottom=0,
            top=np.histogram(arr, bins=bins)[0]/len(arr),
            left=bins[:-1],
            right=bins[1:],
            line_width=0
        )
        row.append(fig)
    rows.append(row)

title_row = []
for nice_col in nice_cols:
    fig = bokeh.plotting.figure(width=1200, height=300, title=nice_col)
    fig.title.text_font_size = '72px'
    fig.outline_line_color = None
    title_row.append(fig)
rows.insert(0, title_row)

fig = bokeh.layouts.gridplot(rows)

bokeh.io.export_png(fig, filename=f'{args.outdir}/locuscompare_hist.png')
bokeh.io.export_svg(fig, filename=f'{args.outdir}/locuscompare_hist.svg')

