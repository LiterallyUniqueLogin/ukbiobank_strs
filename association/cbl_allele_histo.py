#!/usr/bin/env python3

import os

import bokeh.io
import bokeh.models
import bokeh.palettes
import bokeh.plotting
import cyvcf2
import numpy as np
import polars as pl
from statsmodels.stats.weightstats import DescrStatsW

import sample_utils

ukb = os.environ['UKB']

chrom = 11
pos = 119077000
repeat_len = 3
var = next(cyvcf2.VCF(f'{ukb}/str_imputed/runs/first_pass/vcfs/annotated_strs/chr{chrom}.vcf.gz')(f'{chrom}:{pos}'))
lens = sorted(np.unique([round(len(var.REF)/repeat_len)] + [round(len(allele)/repeat_len) for allele in var.ALT]))

imperfect_alleles = pl.read_csv(
    f'{ukb}/association/specific_alleles.tab',
    sep='\t'
).filter(
    (pl.col('pos') == pos) & (pl.col('chrom') == chrom) & (pl.col('name') == 'second C to T')
)['alleles'].to_numpy()[0].split(',')

ethnicities = ('white_brits', 'black', 'south_asian', 'chinese')#, 'irish', 'white_other')
ethnic_gts = {}
for ethnicity in ethnicities:
    ethnic_gts[ethnicity] = {len_: [0,0] for len_ in np.arange(min(lens), max(lens) + 1)}

for ethnicity in ethnicities:
    samp_idx = sample_utils.get_samples_idx_ethnicity(ethnicity)

    len_ = round(len(var.REF)/repeat_len)
    imperfect = 1 if 0 in imperfect_alleles else 0
    ethnic_gts[ethnicity][len_][imperfect] += np.sum((
        np.maximum(0, 1 - np.sum(var.format('AP1'), axis=1)) +
        np.maximum(0, 1 - np.sum(var.format('AP2'), axis=1))
    )[samp_idx])/(2*np.sum(samp_idx))

    for idx, seq in enumerate(var.ALT):
        imperfect = 1 if str(idx + 1) in imperfect_alleles else 0
        len_ = round(len(seq)/repeat_len)
        ethnic_gts[ethnicity][len_][imperfect] += np.sum((var.format('AP1') + var.format('AP2'))[samp_idx, idx])/(2*np.sum(samp_idx))

    corr_data = np.full((len(lens)*2, 2), np.nan)
    corr_data[::2, 0] = lens
    corr_data[1::2, 0] = lens
    corr_data[:, 1] = [0, 1]*len(lens)
    corr_weights = np.full(len(lens)*2, np.nan)
    corr_weights[::2] = [ethnic_gts[ethnicity][len_][0] for len_ in lens]
    corr_weights[1::2] = [ethnic_gts[ethnicity][len_][1] for len_ in lens]
    r2 = DescrStatsW(corr_data, corr_weights).corrcoef[0,1]**2
    print(ethnicity, f'{r2:.3}')
    

max_len = max(len_ for len_ in lens if any(sum(ethnic_gts[ethnicity][len_]) > 0.01 for ethnicity in ethnicities))
min_len = min(len_ for len_ in lens if any(sum(ethnic_gts[ethnicity][len_]) > 0.01 for ethnicity in ethnicities))
cats = [(str(len_), ethnicity) for len_ in range(min_len, max_len + 1) for ethnicity in ethnicities]

fig = bokeh.plotting.figure(
    width= 1600,
    height = 900,
    y_axis_label = 'allele frequency',
    x_axis_label = 'length alleles',
    x_range = bokeh.models.FactorRange(*cats),
    output_backend = 'svg'
)
fig.axis.axis_label_text_font_size = '30px'
fig.axis.major_label_text_font_size = '24px'
fig.xaxis.group_text_font_size = '24px'
fig.xaxis.subgroup_text_font_size = '24px'
fig.xgrid.grid_line_color = None
fig.xaxis.major_label_text_color = None
fig.xaxis.major_tick_line_color = None
fig.toolbar_location = None

cds = bokeh.models.ColumnDataSource(dict(
    x=cats,
    subcats=list(eth.replace('brits', 'british').replace('_', ' ').title() for eth in ethnicities)*(max_len - min_len + 1),
    perfect=[ethnic_gts[ethnicity][int(len_)][0] for (len_, ethnicity) in cats],
    imperfect=[ethnic_gts[ethnicity][int(len_)][1] for (len_, ethnicity) in cats],
))
types = ['perfect', 'imperfect']
fig.vbar_stack(
    types,
    x='x',
    width = 0.9,
    hatch_pattern = [' ', '/'],
    fill_color = bokeh.transform.factor_cmap('x', palette=bokeh.palettes.Colorblind[len(ethnicities)], factors=ethnicities, start=1, end=2),
    legend_group = 'subcats',
    source=cds
)

# remove the legends for  imperfections
for _ in ethnicities:
    del fig.legend[0].items[len(ethnicities)]
fig.legend[0].items.insert(0, fig.legend[0].items[3])
del fig.legend[0].items[4]
fig.legend[0].items.insert(4, fig.legend[0].items[2])
del fig.legend[0].items[2]

# plot again, but legend for imperfections
# since the plots perfectly overlap, shouldn't be an issue
fig.vbar_stack(
    types,
    x='x',
    width = 0.9,
    hatch_pattern = [' ', '/'],
    fill_color = bokeh.transform.factor_cmap('x', palette=bokeh.palettes.Colorblind[len(ethnicities)], factors=ethnicities, start=1, end=2),
    legend_label = types,
    source=cds
)

fig.legend.label_text_font_size = '30px'

bokeh.io.export_png(fig, filename=f'{ukb}/association/allele_plots/cbl.png')
bokeh.io.export_svg(fig, filename=f'{ukb}/association/allele_plots/cbl.svg')
