#!/usr/bin/env python3

import argparse
import ast
import os

import bokeh.embed
import bokeh.events
import bokeh.models
import bokeh.models.callbacks
import bokeh.plotting
import bokeh.resources
import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import WLS
from statsmodels.stats.weightstats import DescrStatsW

import python_array_utils as utils

ukb = os.environ['UKB']

def prep_dict(d):
    return {float(k): v for k, v in d.items() if not np.isnan(v)}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    args = parser.parse_args()
    phenotype = args.phenotype

    title = phenotype.capitalize() + ' signals'

    signals = pd.read_csv(
        f'{ukb}/finemapping/finemap_results/{phenotype}/summary/all_STR_contribs.tab',
        skiprows = 2,
        delimiter='\t'
    )
    signal_split = signals['signal'].str.split('_', n=1, expand=True)
    signals['chrom'] = signal_split[0].astype(int)
    str_split = signals['STR'].str.split('_', n=1, expand=True)
    signals['pos'] = str_split[1].astype(int)

    str_results_fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    associations = pd.read_csv(
        str_results_fname,
        header=0,
        delimiter='\t',
        dtype=utils.get_dtypes(str_results_fname, {'locus_filtered': str})
    )

    signals = signals.merge(
        associations,
        on=['chrom', 'pos'],
        how='inner'
    )
    signals = signals[signals[f'p_{phenotype}'] < 5e-8]
    nrows = signals.shape[0]
    data_dict = {}
    cols = {
        'chrom',
        'pos',
        'pcausal',
        'subset_total_per_allele_dosages',
        f'mean_{phenotype}_per_single_dosage'
    }
    for col in cols:
        data_dict[col] = signals[col].to_numpy()
    del signals

    # calculate effect sizes and third largest allele freq
    data_dict['effect_size'] = np.full(nrows, np.nan)
    data_dict['third_on_allele_freqs'] = np.full(nrows, np.nan)
    for idx in range(nrows):
        single_dosages = prep_dict(ast.literal_eval(data_dict['subset_total_per_allele_dosages'][idx]))
        total_dosage = sum(single_dosages.values())
        single_freqs = {key: value/total_dosage for key, value in single_dosages.items()}

        if len(single_freqs) <= 2:
            data_dict['third_on_allele_freqs'][idx] = 0
        else:
            data_dict['third_on_allele_freqs'][idx] = sum(sorted(single_freqs.values())[:-2])

        doubled_freqs = {}
        for k1 in single_freqs:
            for k2 in single_freqs:
                key = k1 + k2
                if key not in doubled_freqs:
                    doubled_freqs[key] = 0
                doubled_freqs[key] += single_freqs[k1]*single_freqs[k2]

        phenotype_vals = prep_dict(ast.literal_eval(data_dict[f'mean_{phenotype}_per_single_dosage'][idx]))
        for key in phenotype_vals:
            if key not in doubled_freqs:
                print(phenotype_vals)
                print(doubled_freqs)
                assert False
        for key in doubled_freqs:
            if key not in phenotype_vals:
                if not doubled_freqs[key] <= 0.00001:
                    print(phenotype_vals)
                    print(doubled_freqs)
                    assert False

        xs = sorted(phenotype_vals.keys())
        ys = np.array([phenotype_vals[x] for x in xs])
        weights = np.array([doubled_freqs[x] for x in xs])
        xs = np.array(xs).reshape(-1, 1)
        xs = np.concatenate((xs, np.ones((xs.shape[0], 1))), axis=1)

        model = WLS(ys, xs, weights)
        reg_result = model.fit()
        coef = reg_result.params[0]
        scaled_coef = coef*DescrStatsW(xs[:, 0], weights=weights).std
        data_dict['effect_size'][idx] = scaled_coef
    print(data_dict['effect_size'])
    print(data_dict['third_on_allele_freqs'])

    source = bokeh.models.ColumnDataSource(data_dict)

    tooltips = [
        ('chrom', '@chrom'),
        ('pos', '@pos'),
        ('posterior prob causality', '@pcausal'),
        ('multiallelicness', '@third_on_allele_freqs'),
        ('effect size per standard deviation change in repeat size', '@effect_size'),
    ]

    plot = bokeh.plotting.figure(
        width=800,
        height=800,
        title=title,
        x_axis_label='Multiallelicness (1 - frequencies of first and second most common alleles)',
        y_axis_label='Posterior prob causality',
        tools='save,hover',
        tooltips=tooltips
    )
    wheel_zoom = bokeh.models.tools.WheelZoomTool()
    plot.add_tools(wheel_zoom)
    plot.toolbar.active_scroll = wheel_zoom

    pan = bokeh.models.tools.PanTool()
    plot.add_tools(pan)
    plot.toolbar.active_drag = pan

    one_tooltip_callback = bokeh.models.callbacks.CustomJS(code="""
        if (window['one_tooltip']) {
            return;
        }
        window['one_tooltip'] = true;

        var styles = `
            div.bk-tooltip.bk-right>div.bk>div:not(:first-child) {
                display:none !important;
            }
            div.bk-tooltip.bk-left>div.bk>div:not(:first-child) {
                display:none !important;
            }
        `
        var styleSheet = document.createElement("style")
        styleSheet.type = "text/css"
        styleSheet.innerText = styles
        document.head.appendChild(styleSheet)
    """)
    plot.js_on_event(bokeh.events.MouseEnter, one_tooltip_callback)

    display_size = np.abs(data_dict['effect_size'])
    # this puts display sizes dependent on effect_size as an exponential scale between 4 and 10
    display_size -= np.min(display_size)
    display_size *= (np.log(10) - np.log(4))/np.max(display_size)
    display_size += np.log(4)
    display_size = np.exp(display_size)
    source.data['display_size'] = display_size
    plot.circle(
        'third_on_allele_freqs',
        'pcausal',
        source=source,
        size='display_size'
    )

    html = bokeh.embed.file_html(plot, bokeh.resources.CDN, title)
    with open(f'{ukb}/finemapping/finemap_results/{phenotype}/summary/signal_strengths.html', 'w') as outfile:
        outfile.write(html)

if __name__ == '__main__':
        main()
