#!/usr/bin/env python3

import os

import bokeh.io
import bokeh.models
import bokeh.plotting
import numpy as np

import PACSIN2_gts
import sample_utils


ukb = os.environ['UKB']

for name in [
    'dosages',
    'hardcalls',
]:
    sample_idx = sample_utils.get_samples_idx_ethnicity('white_brits')
    idx_gts = np.load(f'{ukb}/finemapping/PACSIN2/idx_{name}.npy')[sample_idx, :]
    allele_details = PACSIN2_gts.allele_details.copy()[:, 1:]
    #allele_details = PACSIN2_gts.allele_details.copy()
    allele_details = np.hstack((allele_details, np.arange(35).reshape(-1, 1)))
    sort_order = np.array(sorted(allele_details.tolist()))[:, -1]
    widths = []
    for allele in range(35):
        if name == 'dosages':
            widths.append(np.sum(idx_gts[:, :, allele]))
        else:
            assert name == 'hardcalls'
            widths.append(np.sum(idx_gts == allele))
    widths = np.array([0] + list(np.array(widths)[sort_order]))

    fig = bokeh.plotting.Figure(
        width=1200,
        height=600,
        title=f'Complex repeat haplotype {name}',
        #y_axis_label='repeat length',
        y_axis_label='number of base pairs (ignoring whitespace)',
        x_axis_label='Relative width implies allele frequency'
    )
    fig.grid.grid_line_color = None
    fig.axis.major_tick_line_color = None
    fig.axis.minor_tick_line_color = None
    fig.axis.major_label_text_color = None
    cummulative_heights = np.zeros(35)
    #for len_, color in zip(range(1, 5), ('#E35669', '#9831E0', '#65D4E6', '#6EF719')):
    for len_, color in zip(range(4), ('#E35669', '#9831E0', '#65D4E6', '#6EF719')):
        varname = PACSIN2_gts.varnames[len_ + 1][8:]
        fig.quad(
            bottom=cummulative_heights[sort_order],
            top=cummulative_heights[sort_order] + allele_details[sort_order, len_],
            left=np.cumsum(widths[:-1]),
            right=np.cumsum(widths[1:]),
            color=[color if
                   (len_ != 2 or snp_num != 1) and
                   (len_ != 1 or snp_num != 2) and
                   (len_ not in {2, 3} or snp_num != 3) and
                   (len_ != 2 or snp_num != 4) and
                   (len_ != 3 or snp_num != 5)
                   else '#FFFF00'
                   for snp_num in allele_details[sort_order, -2]],
            legend_label=varname
        )
        fig.quad(
            bottom=[0],
            top=[0],
            left=[0],
            right=[0],
            color='#FFFF00',
            legend_label='Repeat variant with SNP'
        )
        cummulative_heights += np.max(allele_details[:, len_]) # allele_details[:, len_]
    bokeh.io.export_png(fig, filename=f'{ukb}/finemapping/PACSIN2/white_brits_{name}_Seattle_plot.png')

