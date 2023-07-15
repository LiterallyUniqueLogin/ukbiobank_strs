#!/usr/bin/env python3

import argparse

import bokeh.io
import bokeh.models
import bokeh.palettes
import bokeh.plotting
import bokeh.transform
import cyvcf2
import numpy as np

import trtools.utils.tr_harmonizer as trh

parser = argparse.ArgumentParser()
#parser.add_argument('vcf')
#parser.add_argument('samples', nargs=4)
parser.add_argument('outprefix')
args = parser.parse_args()

#chrom = 'chr14'
#pos = 64247333
#vcf = cyvcf2.VCF(args.vcf)
#
#sampless = []
#for sample_list in args.samples:
#    sampless.append(np.isin(vcf.samples, [line.split('_')[0] for line in open(sample_list).readlines()]))
#    assert np.sum(sampless[-1]) > 100
#
#var = next(vcf(f'{chrom}:{pos}'))
#rec = trh.HarmonizeRecord('hipstr', var)
#gts = rec.GetLengthGenotypes()[:, :-1]

ethnicities = ('white_brits', 'black', 'south_asian', 'chinese')
#ethnic_freqs = {}
#for ethnicity, samples in zip(ethnicities, sampless):
#    ethnic_freqs[ethnicity] = {}
#    ethnic_gts = gts[samples, :]
#    ethnic_gts = ethnic_gts[ethnic_gts > 0]
#    lengths = np.unique(ethnic_gts, return_counts=True)[0]
#    freqs = np.unique(ethnic_gts, return_counts=True)[1]/(ethnic_gts.size)
#    for length, freq in zip(lengths, freqs):
#        ethnic_freqs[ethnicity][length] = freq
#        if length not in [5,6,7]:
#            assert freq < 1e-3


ethnic_freqs = {'white_brits': {5.0: 0.07171309968072935, 6.0: 0.44764380656870906, 7.0: 0.4805801980316635, 7.25: 5.690565043158443e-05, 8.0: 5.990068466482572e-06}, 'black': {5.0: 0.5528428093645485, 5.75: 0.00016722408026755852, 6.0: 0.09882943143812709, 7.0: 0.34816053511705686}, 'south_asian': {4.0: 0.0001497454327643007, 5.0: 0.1305780173704702, 5.75: 0.0002994908655286014, 6.0: 0.5329439952081462, 7.0: 0.3360287511230907}, 'chinese': {5.0: 0.1861788617886179, 6.0: 0.3040650406504065, 7.0: 0.5097560975609756}}

max_len = 7
min_len = 5
cats = [(str(len_), ethnicity) for len_ in range(min_len, max_len + 1) for ethnicity in ethnicities]

fig = bokeh.plotting.figure(
    width= 1600,
    height = 900,
    y_axis_label = 'Allele frequency',
    x_axis_label = 'Allele length',
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
    alleles=[ethnic_freqs[ethnicity][int(len_)] for (len_, ethnicity) in cats],
))
fig.vbar(
    x='x',
    top='alleles',
    width = 0.9,
    fill_color = bokeh.transform.factor_cmap('x', palette=bokeh.palettes.Colorblind[len(ethnicities)], factors=ethnicities, start=1, end=2),
    legend_group = 'subcats',
    source=cds
)

fig.legend.label_text_font_size = '30px'

fig.legend[0].items.insert(0, fig.legend[0].items[3])
del fig.legend[0].items[4]
fig.legend[0].items.insert(2, fig.legend[0].items[3])
del fig.legend[0].items[4]

bokeh.io.export_png(fig, filename=f'{args.outprefix}.png')
bokeh.io.export_svg(fig, filename=f'{args.outprefix}.svg')

