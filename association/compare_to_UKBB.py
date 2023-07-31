#!/usr/bin/env python3

import argparse
import gzip

import bokeh.io
import bokeh.models
import bokeh.models.formatters
import bokeh.models.tickers
import bokeh.plotting
import bokeh.transform
import matplotlib.cm
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
#from qqman import qqman

#def validate_our_code():
#    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))
#
#    print('loading csvs ... ', flush=True)
#    my_df = pd.read_csv(
#        f'{ukb}/association/plots/input/eosinophil_count/my_imputed_snp_chr21_results.tab',
#        header=0,
#        delimiter='\t',
#        usecols=['chrom', 'pos', 'alleles', 'p_eosinophil_count', 'locus_filtered'],
#        dtype={'chrom': int, 'pos': int, 'alleles': str, 'p_eosinophil_count': float, 'locus_filtered': str},
#    )
#    my_df.rename(columns={'chrom': 'CHR', 'pos': 'BP', 'p_eosinophil_count': 'P'}, inplace=True)
#    alleles = my_df['alleles'].str.split(',', 1, expand=True)
#    print(my_df.shape)
#    
#    my_df = my_df[my_df['locus_filtered'] == 'False']
#    print(my_df.shape)
#    my_df['REF'] = alleles.iloc[:, 0]
#    my_df['ALT'] = alleles.iloc[:, 1]
#
#    plink_df = pd.read_csv(
#        f'{ukb}/association/results/eosinophil_count/plink_snp/results.tab',
#        header=0,
#        delimiter='\t',
#        usecols=['#CHROM', 'POS', 'P', 'REF', 'ALT'],
#        dtype={'#CHROM': int, 'POS': int, 'P': float}
#    )
#    plink_df.rename(columns={'#CHROM': 'CHR', 'POS': 'BP'}, inplace=True)
#    plink_df = plink_df[plink_df['CHR'] == 21]
#    print(plink_df.shape)
#
#    plink_df = plink_df.merge(
#        my_df,
#        how='inner',
#        on=['BP', 'REF', 'ALT'],
#        suffixes=(None, '_mine')
#    )
#    plink_df = plink_df[['CHR', 'BP', 'P', 'REF', 'ALT']]
#    print(plink_df.shape)
#
#    my_df = my_df.merge(
#        plink_df,
#        how='inner',
#        on=['BP', 'REF', 'ALT'],
#        suffixes=(None, '_theirs')
#    )
#    my_df = my_df[['CHR', 'BP', 'P']]
#    print(my_df.shape)
#
#    print('plotting manhattans ... ', flush=True)
#    qqman.manhattan(
#        my_df,
#        cmap=matplotlib.cm.get_cmap('viridis'),
#        ax=ax1,
#        title="P-values from this study's code",
#        suggestiveline=False
#    )
#
#    qqman.manhattan(
#        plink_df,
#        cmap=matplotlib.cm.get_cmap('viridis'),
#        ax=ax2,
#        title='Plink p-values',
#        suggestiveline=False
#    )
#
#    plt.savefig(f'{ukb}/export_scripts/results/validate_our_code.png')

def linear_int_interpolate(c1, c2, dist):
    c_new = []
    for coord1, coord2 in zip(c1, c2):
        c_new.append(coord1 + round((coord2 - coord1)*dist))
    return c_new

def scatter_with_panukbb(outdir, pan_ukbb_tsv, my_pipeline_tsv):
    print('loading panukbb ... ', flush=True)
    panukbb_df = pl.scan_csv(
        #f'{ukb}/misc_data/snp_summary_stats/bilirubin/neale/biomarkers-30840-both_sexes-irnt.tsv',
        pan_ukbb_tsv,
        sep='\t',
        dtypes={'chr': str, 'pos': int, 'ref': str, 'alt': str, 'pval_EUR': float},
        null_values='NA'
    ).rename({'pval_EUR': 'p'}).filter(
        pl.col('chr') != 'X'
    ).with_column(
        pl.col('chr').cast(int)
    ).select(['chr', 'pos', 'ref', 'alt', 'p'])

    print('loading plink ... ', flush=True)
    my_pipeline_df = pl.scan_csv(
        my_pipeline_tsv,
        sep='\t',
        dtypes={'#CHROM': int, 'POS': int, 'REF': str, 'ALT': str, 'P': float},
        null_values='NA'
    ).rename(
        {'#CHROM': 'chr', 'POS': 'pos', 'REF': 'ref', 'ALT': 'alt', 'P': 'p'}
    ).select(['chr', 'pos', 'ref', 'alt', 'p'])

    print('merging ... ', flush=True)
    n_rects = 100
    df = panukbb_df.join(
        my_pipeline_df,
        how='inner',
        on=['chr', 'pos', 'ref', 'alt'],
        suffix='_mine'
    ).with_columns([
        (-pl.max([1e-50, pl.col('p')]).log10()).alias('p'),
        (-pl.max([1e-50, pl.col('p_mine')]).log10()).alias('p_mine'),
    ]).select([
        (pl.col('p')*n_rects).floor()/n_rects,
        (pl.col('p_mine')*n_rects).floor()/n_rects
    ]).groupby(['p', 'p_mine']).agg(pl.count()).collect()

    print('plotting ... ', flush=True)

    fig_height=1200
    fig = bokeh.plotting.figure(
        width=fig_height,
        height=fig_height,
        y_axis_label='-log10(p-value our pipeline)',
        x_axis_label='-log10(p-value Pan UKBB)',
        x_range=[0,50],
        y_range=[0,50],
    )
    fig.axis.axis_label_text_font_size = '36px'
    fig.axis.major_label_text_font_size = '30px'
    palette = [linear_int_interpolate((134,204,195), (9,41,46), i/254) for i in range(-1, 255)]
    cmap = bokeh.transform.log_cmap(
        'count',
        palette=palette,
        low=1,
        high=max(df['count'].to_numpy()),
        low_color=(255, 255, 255)
    )
    color_mapper = bokeh.models.LogColorMapper(
        palette = palette,
        low=0,
        high=max(df['count'].to_numpy())
    )

    cds = bokeh.models.ColumnDataSource(dict(
        left=df['p'].to_numpy(), right=df['p'].to_numpy() + 50/n_rects,
        bottom=df['p_mine'].to_numpy(), top=df['p_mine'].to_numpy() + 50/n_rects,
        count=df['count'].to_numpy()
    ))
    fig.quad(
        left='left', right='right', bottom='bottom', top='top', source=cds, fill_color=cmap, line_width=0
    )
    fig.line(
        x=[0,50],
        y=[0,50],
        line_width=4,
        color='black',
        line_dash = 'dashed'
    )

    color_bar = bokeh.models.ColorBar(
        title = 'Number of loci',
        color_mapper = cmap['transform'],
        width=35,
        height=fig_height//3,
        title_text_font_size = '30px',
        major_label_text_font_size = '30px',
        ticker = bokeh.models.tickers.FixedTicker(ticks=[0, 10, 100, 1000, 10000, 100000]),
        formatter = bokeh.models.formatters.NumeralTickFormatter(format = '0a'),
        location='bottom_right'
    )
    fig.add_layout(color_bar, 'right')
    fig.toolbar_location = None
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color=None
    bokeh.io.export_svg(fig, filename=f'{outdir}/panukbb_scatter.svg')
    bokeh.io.export_png(fig, filename=f'{outdir}/panukbb_scatter.png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('pan_ukbb_tsv')
    parser.add_argument('my_pipeline_tsv')
    args = parser.parse_args()
    #validate_our_code()
    scatter_with_panukbb(args.outdir, args.pan_ukbb_tsv, args.my_pipeline_tsv)
