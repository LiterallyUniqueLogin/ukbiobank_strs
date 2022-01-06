#!/usr/bin/env python3

import argparse
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

    title = phenotype.capitalize() + ' signals (size measures signal strength)'

    signals = pd.read_csv(
        f'{ukb}/finemapping/summary/{phenotype}_table.tab',
        delimiter='\t',
        encoding='UTF-8',
        dtype={'included_only_due_to_literature': str, 'included_only_due_to_curation': str}
    )
    signals = signals[(signals['included_only_due_to_literature'] != 'True') & (signals['included_only_due_to_curation'] != 'True')]
    data_dict = {}
    cols = {
        'chrom',
        'start_pos',
        'FINEMAP_pcausal',
        'subset_multiallelicness'
    }
    for col in cols:
        data_dict[col] = signals[col].to_numpy()
    data_dict['effect_size'] = signals['Δphenotype_per_s.d._increase_in_repeat_size'].to_numpy()

    del signals

    source = bokeh.models.ColumnDataSource(data_dict)
    print(source)

    tooltips = [
        ('chrom', '@chrom'),
        ('pos', '@start_pos'),
        ('posterior prob causality', '@FINEMAP_pcausal'),
        ('multiallelicness', '@subset_multiallelicness'),
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
    # this puts display sizes dependent on effect_size as an exponential scale between 4 and 20
    display_size -= np.min(display_size)
    display_size *= (np.log(20) - np.log(4))/np.max(display_size)
    display_size += np.log(4)
    display_size = np.exp(display_size)
    source.data['display_size'] = display_size
    plot.circle(
        'subset_multiallelicness',
        'FINEMAP_pcausal',
        source=source,
        size='display_size',
        alpha=0.8
    )

    html = bokeh.embed.file_html(plot, bokeh.resources.CDN, title)
    with open(f'{ukb}/finemapping/finemap_results/{phenotype}/summary/signal_strengths.html', 'w') as outfile:
        outfile.write(html)

if __name__ == '__main__':
        main()
