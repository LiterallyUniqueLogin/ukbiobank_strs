#!/usr/bin/env python3

import argparse

import bokeh.io
import bokeh.models
import bokeh.models.tickers
import bokeh.plotting
import numpy as np

other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']

parser = argparse.ArgumentParser()
parser.add_argument('outdir')
parser.add_argument('enrichment_stats')
args = parser.parse_args()

def barplot_fig(cats, data, ps, legend):
    cats = [(
        cat[0].replace('_', ' ').
               replace(' trinucs', '\ntrinucleotides').
               replace('non ', 'non-').
               replace('transcribed ', 'transcribed\n').
               replace('upstream ', 'upstream\n').
               replace(' dinucs', '\ndinucleotides').
               replace('tetranucs', 'tetra-\nnucleotides').
               replace('pentanucs', 'penta-\nnucleotides').
               replace('hexanucs', 'hexa-\nnucleotides').
               replace('unclear ', 'unclear\n'),
        cat[1]
    ) for cat in cats]
    fig = bokeh.plotting.figure(
        y_axis_label = 'Prevalence (%)',
        x_range = bokeh.models.FactorRange(*cats),
        width=len(data)//3*200,
        height=800,
        toolbar_location=None,
        output_backend='svg'
    )
    fig.xaxis.major_label_text_color = None
    fig.xaxis.major_tick_line_color = None
    fig.xaxis.group_text_font_size = '24px'
    fig.xaxis.subgroup_text_font_size = '24px'
    fig.axis.axis_label_text_font_size = '36px'
    fig.axis.major_label_text_font_size = '30px'
    fig.title.text_font_size = '36px'
    arr_data = np.array(data)
    arr_data[arr_data == 0] = np.nan
    condition_names = np.unique(list(cat[1] for cat in cats))
    cds = bokeh.models.ColumnDataSource(dict(
        x=cats,
        top=arr_data,
        subcats=[cat[1] for cat in cats]
    ))
    kwargs = dict(
        x='x',
        top='top',
        source=cds,
        width=0.9,
        fill_color = bokeh.transform.factor_cmap('x', palette=bokeh.palettes.Colorblind[3], factors=condition_names, start=1, end=2),
    )
    if legend:
        kwargs['legend_group'] = 'subcats'
    fig.vbar(
        **kwargs
    )
    if legend:
        fig.legend.label_text_font_size = '30px'
        fig.legend.location = 'top_left'
        fig.legend[0].items.insert(1, fig.legend[0].items[2])
        del fig.legend[0].items[3]
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.x_range.range_padding = 0.1
    #fig.xaxis.major_label_orientation = 1.3
    fig.xgrid.grid_line_color = None
    y_range = max(data) # since plots are relative to zero
    for group in range(len(cats)//3):
        compare_bars(fig, y_range, cats[group*3], cats[group*3+1], data[group*3], data[group*3+1], data[group*3:group*3+2], ps[group*3], 0, -1, 1)
        compare_bars(fig, y_range, cats[group*3+1], cats[group*3+2], data[group*3+1], data[group*3+2], data[group*3+1:group*3+3], ps[group*3+2], 1, 0, 1)
        compare_bars(fig, y_range, cats[group*3], cats[group*3+2], data[group*3], data[group*3+2], data[group*3:group*3+3], ps[group*3+1], -1, 1, 2)
    return fig

def compare_bars(fig, y_range, x1, x2, y1, y2, ys, p_val, x1step, x2step, ystep):
    if p_val > 0.05:
        return
    x1 = [*x1, .15*x1step]
    x2 = [*x2, .15*x2step]
    y_step_size = 0.05*y_range
    top = max(ys) + ystep*y_step_size
    fig.line(x=[x1, x1], y=[y1, top])
    fig.line(x=[x2, x2], y=[y2, top])
    fig.line(x=[x1, x2], y=[top, top])
    if p_val >= 1e-300:
        text=f'p={p_val:.1g}'
    else:
        text='p<1e-300'
    if p_val < 0.05/50:
        text = text + '**'
    elif p_val < 0.05:
        text = text + '*'

    cds = bokeh.models.ColumnDataSource(dict(
        x=[x1[:2]],
        y=[top+0.01*y_step_size],
        text=[text]
    ))
    fig.add_layout(bokeh.models.LabelSet(x='x', y='y', x_offset=-20, text='text', source=cds, text_font_size='20px'))

datas = {}
for fig_loc in 'regions', 'repeats':
    datas[fig_loc] = {}
    for bar_loc in 'left', 'right':
        datas[fig_loc][bar_loc] = {}
        datas[fig_loc][bar_loc]['data'] = []
        datas[fig_loc][bar_loc]['cat'] = []
        datas[fig_loc][bar_loc]['ps'] = []

#f'{ukb}/post_finemapping/results/enrichments.tab') as stats:
with open(args.enrichment_stats) as stats:
    next(stats) #skip header
    for line in stats:
        if 'CCG' in line:
            continue
        line = line.replace('\\n', '\n')
        split = line.split('\t')
        category = split[0]
        ori_perc = float(split[1].split()[1][1:-2])
        if 'other' in category or 'poly' in category or 'unclear' in category:
            fig_loc = 'repeats'
        else:
            fig_loc = 'regions'
        if ori_perc >= 10:
            bar_loc = 'right'
        else:
            bar_loc = 'left'
        curr_data = datas[fig_loc][bar_loc]['data']
        curr_cat = datas[fig_loc][bar_loc]['cat']
        curr_ps = datas[fig_loc][bar_loc]['ps']
        curr_cat.extend([
            (category, f'all STRs'),
            (category, f'genome-wide significant STRs'),
            (category, 'confidently fine-mapped STRs')
        ])
        curr_ps.extend([float(split[idx]) for idx in range(4, 7)])

        curr_data.append(ori_perc)
        curr_data.append(float(split[2].split()[1][1:-2]))
        curr_data.append(float(split[3].split()[1][1:-2]))

fig = bokeh.layouts.column([
    barplot_fig(datas['regions']['left']['cat'], datas['regions']['left']['data'], datas['regions']['left']['ps'], True),
    barplot_fig(datas['regions']['right']['cat'], datas['regions']['right']['data'], datas['regions']['right']['ps'], False),
    barplot_fig(datas['repeats']['left']['cat'], datas['repeats']['left']['data'], datas['repeats']['left']['ps'], True),
    barplot_fig(datas['repeats']['right']['cat'], datas['repeats']['right']['data'], datas['repeats']['right']['ps'], False),
])
bokeh.io.export_svg(fig, filename=f'{args.outdir}/enrichment_barplots.svg')
bokeh.io.export_png(fig, filename=f'{args.outdir}/enrichment_barplots.png')

