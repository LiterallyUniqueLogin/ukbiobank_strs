#!/usr/bin/env python3

import os

import bokeh.io
import bokeh.models
import bokeh.models.formatters
import bokeh.models.tickers
import bokeh.plotting
import bokeh.transform

import numpy as np
import polars as pl

import matplotlib.pyplot as plt
import matplotlib.colors as colors

ukb = os.environ['UKB']

def loadFreqs(freqLocation, identifyingVCFLocation, cNum, vcfHeaderLines):
    freqs = pl.scan_csv(
        f'{freqLocation}/chr{cNum}.afreq',
        sep='\t',
        dtypes={'ALT_FREQS': str}
    ).collect()
    vcf = pl.scan_csv(
        f'{identifyingVCFLocation}/chr{cNum}.vcf',
        skip_rows = vcfHeaderLines,
        sep='\t',
    ).select(
        ['#CHROM', 'POS', 'ID', 'REF', 'ALT']
    ).collect()
    assert np.all(vcf['ID'].to_numpy() == freqs['ID'].to_numpy())
    return vcf.with_column(
        freqs['ALT_FREQS'].alias('afreq')
    ).drop('ID')


def linear_int_interpolate(c1, c2, dist):
    c_new = []
    for coord1, coord2 in zip(c1, c2):
        c_new.append(coord1 + round((coord2 - coord1)*dist))
    return c_new

def graph_entire(df):
    fig_height=1200
    fig = bokeh.plotting.figure(
        width=fig_height,
        height=fig_height,
        y_axis_label='UKB allele frequencies',
        x_axis_label='STR panel allele frequencies',
        x_range=[0,1],
        y_range=[0,1],
        output_backend='svg',
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
        left=df['afreq_snpstr'].to_numpy(), right=df['afreq_snpstr'].to_numpy() + 1/n_rects,
        bottom=df['afreq'].to_numpy(), top=df['afreq'].to_numpy() + 1/n_rects,
        count=df['count'].to_numpy()
    ))
    fig.quad(
        left='left', right='right', bottom='bottom', top='top', source=cds, fill_color=cmap, line_width=0
    )
    fig.line(
        x=[0,1-.12],
        y=[.12,1],
        line_width=4,
        color='red',
        line_dash = 'dashed'
    )
    fig.line(
        x=[.12,1],
        y=[0,1-.12],
        line_width=4,
        color='red',
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
    bokeh.io.export_svg(fig, filename=f'{ukb}/pre_imputation_qc/allele_freqs/contrast.svg')
    bokeh.io.export_png(fig, filename=f'{ukb}/pre_imputation_qc/allele_freqs/contrast.png')

    '''
    graph_xs = xs[:variantCount]
    graph_ys = ys[:variantCount]
    
    line = np.linspace(0, 1, 1000)
    axes.plot(line, line)
    axes.plot(line, a*line**2 + b*line + c, 'r', label="Proposed filtering cutoff")
    axes.plot(line, a2*line**2 + b2*line + c2, 'r')
    axes.scatter(*zip(*filteredSNPs.values()), c='y')
    fig.colorbar(hist, label='count in bin')
    axes.legend()
    axes.set_title(title)
    plt.show()
'''


ukb_freqs = []
snpstr_freqs = []

for cNum in range(1,23):
    print(cNum)
    ukb_freqs.append(loadFreqs(
        f"{ukb}/pre_imputation_qc/allele_freqs/hap_white_brit_high_quality_sane_kinship_unrelated",
        f"{ukb}/microarray/vcf_1_sample",
        cNum, 5
    ))
    snpstr_freqs.append(loadFreqs(
        f"{ukb}/pre_imputation_qc/allele_freqs/snpstr_eur",
        f"{ukb}/snpstr/vcf_1_sample",
        cNum,
        10
    ))
ukb_freqs = pl.concat(ukb_freqs)
snpstr_freqs = pl.concat(snpstr_freqs)

n_rects = 100
merged = ukb_freqs.join(
    snpstr_freqs,
    how = 'inner',
    on = ['#CHROM', 'POS', 'REF', 'ALT'],
    suffix = '_snpstr'
).select([
    (pl.col('afreq').cast(float)*n_rects).floor()/n_rects,
    (pl.col('afreq_snpstr').cast(float)*n_rects).floor()/n_rects
]).groupby(['afreq', 'afreq_snpstr']).agg(pl.count())

'''
print(merged.filter(
    (pl.col('afreq_snpstr').cast(float) - pl.col('afreq').cast(float)).abs() >= .12
).to_pandas())
'''

graph_entire(merged)#'Allele frequency among shared SNPs in post-filtering white European samples', 0, 1, -0.12, 0, 1, 0.12, merged)
