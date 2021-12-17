#!/usr/bin/env python3

import argparse
import os

import bokeh.io
import bokeh.models
import bokeh.models.labeling
import bokeh.models.tickers
import bokeh.models.tools
import bokeh.plotting
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import scipy.stats
import scipy.stats.contingency

import region_plots

ukb = os.environ['UKB']

'''
ALL STRs
And with SIG <= 5e-8
And with FINEMAP pcausal >= .9
'''

parser = argparse.ArgumentParser()
parser.add_argument('ext')
args = parser.parse_args()
ext = args.ext

def graph_relative_rate_by_pcausal(fig, compare_STRs, in_category, category_string, color, dash, ext):
    xs = np.arange(0.00, 1.01, 0.01)
    finemaped_subset_counts = []
    comparison_count = compare_STRs.shape[0]
    comparison_subset_count = pl.sum(in_category)
    rate_ratios = [np.nan]
    upper_cis = [np.nan]
    lower_cis = [np.nan]

    for x in xs:
        finemaped_count = pl.sum(compare_STRs['FINEMAP_pcausal'] >= x)
        finemaped_subset_count = pl.sum(in_category & (compare_STRs['FINEMAP_pcausal'] >= x))
        finemaped_subset_counts.append(finemaped_subset_count)
        if x > 0:
            result = scipy.stats.contingency.relative_risk(
                finemaped_subset_count,
                finemaped_count,
                comparison_subset_count - finemaped_subset_count,
                comparison_count - finemaped_count
            )
            rate_ratios.append(result.relative_risk)
            # +/-1 SD
            lower, upper = result.confidence_interval(.68)
            lower_cis.append(lower)
            upper_cis.append(upper)
    rate_ratios = np.array(rate_ratios)
    finemaped_subset_counts = np.array(finemaped_subset_counts)

    source = bokeh.models.ColumnDataSource(dict(
        xs=xs,
        rate_ratios=rate_ratios
    ))
    line = fig.line(x='xs', y='rate_ratios', source=source, legend_label=category_string, color=color, line_dash=dash, line_width=5)
    if ext == 'html':
        ci = fig.varea(xs, lower_cis, upper_cis, color=color, alpha=0.2)
        ci.visible = False
        return (source, ci)
    return source

compare_STRs = pl.scan_csv(
    f'{ukb}/side_analyses/str_annotations/all_STRs.tab',
    sep='\t'
).filter(pl.col('gwas_sig') & ~pl.col('FINEMAP_pcausal').is_null()).with_column(
    (pl.col('FINEMAP_pcausal') >= 0.9).alias('FINEMAPed')
).collect()

def setup_relative_fig():
    relative_enrichment_fig = bokeh.plotting.figure(
        title='Relative rate of STR categories by FINEMAP posterior causality',
        y_axis_label = 'Relative rate',
        x_axis_label = 'FINEMAP posterior causality threshold',
        width=1200,
        height=800
    )
    relative_enrichment_fig.background_fill_color = None
    relative_enrichment_fig.border_fill_color = None
    relative_enrichment_fig.grid.grid_line_color = None
    relative_enrichment_fig.toolbar_location = None
    relative_enrichment_fig.title.text_font_size = '30px'
    relative_enrichment_fig.axis.axis_label_text_font_size = '26px'
    relative_enrichment_fig.axis.major_label_text_font_size = '20px'
    relative_enrichment_fig.line([-0,1], [1,1], line_dash='solid', line_width=2, color='black')
    ticks = np.arange(0, 1.01, .2)
    relative_enrichment_fig.xaxis.ticker = ticks
    finemaped_counts = [pl.sum(compare_STRs['FINEMAP_pcausal'] >= tick) for tick in ticks]
    relative_enrichment_fig.xaxis.major_label_overrides = {
        f'{tick:.1g}': f'{tick:.1g}\n(n = {count})' for tick, count in zip(ticks, finemaped_counts)
    }
    return relative_enrichment_fig

relative_enrichment_fig = setup_relative_fig()
'''
relative_enrichment_fig.add_layout(bokeh.models.Title(
    text='Relative rate is compared to the set of all non-MHC gwas-sig STRs' #: {comparison_subset_count}/{comparison_count} = {comparison_subset_count/comparison_count*100:.4f}%'
), 'below')
relative_enrichment_fig.add_layout(bokeh.models.Title(
    text='Relative rate is defined as (% category STRs among FINEMAP STRs with posterior causality >= x) / (% category STRs among all non-MHC gwas-sig STRs)'
), 'below')
relative_enrichment_fig.add_layout(bokeh.models.Title(
    text='gwas-sig STRs are those with association p <= 5e-8, not in the MHC. FINEMAP STRs must also be gwas-sig'
), 'below')
'''

plots = []
dashes=('solid', 'dashed', 'dotted')
for category, dash in zip(('exonic', 'intronic', 'transcribed_non_protein'), dashes):
    plots.append(graph_relative_rate_by_pcausal(
        relative_enrichment_fig,
        compare_STRs,
        compare_STRs[category],
        category,
        color='#59B65E',
        dash=dash,
        ext=ext
    ))

plots.append(graph_relative_rate_by_pcausal(
    relative_enrichment_fig,
    compare_STRs,
    compare_STRs['unit'] == 'A',
    'polyA',
    color='#971FED',
    dash='solid',
    ext=ext
))

plots.append(graph_relative_rate_by_pcausal(
    relative_enrichment_fig,
    compare_STRs,
    compare_STRs['unit'] == 'AC',
    'AC dinucleotide',
    color='#971FED',
    dash='dashed',
    ext=ext
))
'''
plots.append(graph_relative_rate_by_pcausal(
    relative_enrichment_fig,
    compare_STRs,
    compare_STRs['period'] == 3,
    'trinucleotide',
    color='#971FED',
    dash='dotted',
    ext=ext
))

plots.append(graph_relative_rate_by_pcausal(
    relative_enrichment_fig,
    compare_STRs,
    (compare_STRs['period'] == '3') & (compare_STRs['exonic']),
    'exonic and trinuc',
    color='#E61A0B',
    dash='dotted',
    ext=ext
))

plots.append(graph_relative_rate_by_pcausal(
    relative_enrichment_fig,
    compare_STRs,
    (compare_STRs['period'] != '3') & (compare_STRs['exonic']),
    'exonic not trinuc',
    color='#E61A0B',
    dash='dashed',
    ext=ext
))
'''

relative_enrichment_fig.legend.label_text_font_size = '22px'
relative_enrichment_fig.legend.location = 'top_left'

if ext == 'html':
    tap = bokeh.models.tools.TapTool()
    relative_enrichment_fig.add_tools(tap)
    relative_enrichment_fig.toolbar.active_tap = tap
    callback = bokeh.models.CustomJS(
        args=dict(
            sources=[plot[0] for plot in plots],
            cis=[plot[1] for plot in plots]
        ),
        code="""
            for(var i=0;i<sources.length;i++) {
                console.log(sources[i].selected.line_indices)
                if (sources[i].selected.line_indices.length %2 == 1) {
                    cis[i].visible=true;
                } else {
                    cis[i].visible=false;
                }
                sources[i].selected.line_indices = [];
            }
        """
    )
    tap.callback = callback

region_plots.export(relative_enrichment_fig, f'{ukb}/post_finemapping/results/FINEMAP_relative_rate.{ext}', ext)

trinuc_exonic_fig = setup_relative_fig()
plots = []
plots.append(graph_relative_rate_by_pcausal(
    trinuc_exonic_fig,
    compare_STRs,
    (compare_STRs['period'] == '3') & (compare_STRs['exonic']),
    'exonic trinucleotide',
    color='#E61A0B',
    dash='solid',
    ext=ext
))

trinuc_exonic_fig.legend.label_text_font_size = '22px'
trinuc_exonic_fig.legend.location = 'top_left'

if ext == 'html':
    tap = bokeh.models.tools.TapTool()
    trinuc_exonic_fig.add_tools(tap)
    trinuc_exonic_fig.toolbar.active_tap = tap
    callback = bokeh.models.CustomJS(
        args=dict(
            sources=[plot[0] for plot in plots],
            cis=[plot[1] for plot in plots]
        ),
        code="""
            for(var i=0;i<sources.length;i++) {
                console.log(sources[i].selected.line_indices)
                if (sources[i].selected.line_indices.length %2 == 1) {
                    cis[i].visible=true;
                } else {
                    cis[i].visible=false;
                }
                sources[i].selected.line_indices = [];
            }
        """
    )
    tap.callback = callback

region_plots.export(trinuc_exonic_fig, f'{ukb}/post_finemapping/results/FINEMAP_trinuc_exonic_relative_rate.{ext}', ext)


'''
fig, ax = plt.subplots()
ax.set_xlim(0, np.max(compare_STRs['mean_len']))
ax.set_title('CDF of STR mean len')
ax.set_ylabel('density')
ax.set_xlabel('distance (bp)')
_BetterCDF(compare_STRs['mean_len'], ax)
p_val = scipy.stats.mannwhitneyu(compare_STRs['mean_len'], compare_STRs.loc[compare_STRs['FINEMAPed'], 'mean_len'])[1]
_BetterCDF(compare_STRs.loc[compare_STRs['FINEMAPed'], 'mean_len'], ax)
legends = ['gwas sig STRs', 'FINEMAPed']
ax.text(0.5, -0.1, "Comparing FINEMAPed to gwas_sig. gwas_sig is association p<=5e-8. FINEMAPed is that and also FINEMAP posterior causal >= 0.9. Both groups exclude MHC", ha="center", transform=ax.transAxes)
ax.text(0.5, -1.1, f"Medians: {np.median(compare_STRs['mean_len'])} (gwas_sig) {np.median(compare_STRs.loc[compare_STRs['FINEMAPed'], 'mean_len'])} (FINEMAPed)", ha="center", transform=ax.transAxes)
plt.legend(legends)
plt.savefig(f'{ukb}/post_finemapping/results/mean_len_cdf.png')
plt.savefig(f'{ukb}/post_finemapping/results/mean_len_cdf.pdf')

for category, class_ in ('gene', 'intergenic'), ('exon', 'intronic'):
    for stream in 'upstream', 'downstream':
        col = f'{stream}_{category}_dist'
        subbed_data = compare_STRs.loc[compare_STRs[class_], :]
        fig, ax = plt.subplots()
        ax.set_xlim(0, np.max(subbed_data[col]))
        ax.set_title(f'CDF of {col} from within {class_} STRs')
        ax.set_ylabel('density')
        ax.set_xlabel('distance (bp)')
        _BetterCDF(subbed_data[col], ax)
        p_val = scipy.stats.mannwhitneyu(subbed_data[col], subbed_data.loc[subbed_data['FINEMAPed'], col])[1]
        _BetterCDF(subbed_data.loc[subbed_data['FINEMAPed'], col], ax)
        legends = ['gwas sig STRs', 'FINEMAPed']
        ax.text(0.5, -0.1, "Comparing FINEMAPed to gwas_sig. gwas_sig is association p<=5e-8. FINEMAPed is that and also FINEMAP posterior causal >= 0.9", ha="center", transform=ax.transAxes)
        ax.text(0.5, -1.1, f"Medians: {np.median(subbed_data[col])} (gwas_sig) {np.median(subbed_data.loc[subbed_data['FINEMAPed'], col])} (FINEMAPed)", ha="center", transform=ax.transAxes)
        plt.legend(legends)
        plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf.png')
        plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf.pdf')
        ax.set_xlim(0, 50000)
        subbed_data = subbed_data.loc[subbed_data[col] < 50000, :]
        p_val = scipy.stats.mannwhitneyu(subbed_data[col], subbed_data.loc[subbed_data['FINEMAPed'], col])[1]
        plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf_50kbp.png')
        plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf_50kbp.pdf')

'''
