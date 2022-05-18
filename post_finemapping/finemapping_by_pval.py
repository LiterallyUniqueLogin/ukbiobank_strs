#!/usr/bin/env python3

import os

import bokeh.io
import bokeh.models
import bokeh.plotting
import numpy as np
import polars as pl

ukb = os.environ['UKB']

df = pl.read_csv(
    f'{ukb}/post_finemapping/results/singly_finemapped_strs_for_paper.tab',
    sep='\t'
).with_column(
    pl.when(
        pl.col('association_p_value') < 1e-300
    ).then(
        pl.lit(300)
    ).otherwise(
        -pl.col('association_p_value').log10()
    ).alias('association_p_value')
)

# Min abs corr across all CSes
fig = bokeh.plotting.figure(
    width=1200,
    height=1200,
    title='Fine-mapped loci by association p-value',
    x_axis_label='-log10(p-value)',
    y_axis_label='Number of fine-mapped loci',
)
fig.axis.axis_label_text_font_size = '30px'
fig.title.text_font_size = '30px'
fig.axis.major_label_text_font_size = '20px'
fig.background_fill_color = None
fig.border_fill_color = None
fig.grid.grid_line_color = None
fig.toolbar_location = None
step = 5
left_edges = np.arange(0, 300 + step, step)
finemap_vals = df.filter(pl.col('finemap_CP') >= 0.8)['association_p_value'].to_numpy().flatten()
finemap_ys = [np.sum((left_edge <= finemap_vals) & (finemap_vals < left_edge + step)) for left_edge in left_edges]
fig.quad(top=finemap_ys, bottom=0, left=left_edges, right=left_edges+step, alpha=0.5, color='blue', legend_label='finemap')
susie_vals = df.filter(pl.col('susie_CP') >= 0.8)['association_p_value'].to_numpy().flatten()
susie_ys = [np.sum((left_edge <= susie_vals) & (susie_vals < left_edge + step)) for left_edge in left_edges]
fig.quad(top=susie_ys, bottom=0, left=left_edges, right=left_edges+step, alpha=0.5, color='red', legend_label='susie')

bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/finemapping_by_pval.png')
bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/finemapping_by_pval.svg')

