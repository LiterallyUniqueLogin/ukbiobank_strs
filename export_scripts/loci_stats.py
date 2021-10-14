#!/usr/bin/env python3

import os

import bokeh.io
import bokeh.plotting
import numpy as np
import pandas as pd
import scipy.stats

ukb = os.environ['UKB']

def main():
    dfs = []
    for chrom in range(1, 23):
        print(f'loading csv {chrom} ...', flush=True)
        dfs.append(pd.read_csv(
            f'{ukb}/export_scripts/intermediate_results/chr{chrom}_loci_summary.tab',
            header=0,
            delimiter='\t',
            dtype={'chr': int, 'pos': int, 'allele_dist': str, 'entropy': float, 'heterozygosity': float, 'multiallelicness': float}
        ))
    df = pd.concat(dfs)

    for col in 'entropy', 'heterozygosity', 'multiallelicness':
        print(f'plotting {col} ...', flush=True)
        max_val = np.max(df[col])
        xs = np.arange(0, max_val*1.001, max_val/1001)
        ys = scipy.stats.gaussian_kde(df[col])(xs)
        x_axis_label = col
        if col == 'entropy':
            x_axis_label += ' (bits)'
        figure = bokeh.plotting.figure(
            title=col.capitalize() + ' distribution',
            x_axis_label=x_axis_label,
            y_axis_label='pdf',
            width=600,
            height=600,
            tools=''
        )
        figure.background_fill_color = None
        figure.border_fill_color = None
        figure.ygrid.grid_line_color = None
        figure.xgrid.grid_line_color = None
        figure.varea(
            x=xs,
            y1=np.zeros(xs.shape),
            y2=ys
        )
        line_x = None
        if col == 'entropy':
            line_x = 1
        elif col == 'heterozygosity':
            line_x = .5
        if line_x:
            figure.line(
                [line_x, line_x], [0, max(ys)], line_dash = 'dotdash', color='red',
                legend_label=f'max biallelic variant {col}'
            )
            
        bokeh.io.export_png(figure, filename=f'{ukb}/export_scripts/results/{col}_distribution.png')
        bokeh.io.export_svg(figure, filename=f'{ukb}/export_scripts/results/{col}_distribution.svg')

    for thresh in (0.001, 0.004, 0.015, 0.06):
        print(f'plotting allele count > {thresh} ... ', flush=True)
        counts = {}
        for str_ in df['allele_dist']:
            split = str_.split(' ')
            count = 0
            for freq in split[1::2]:
                if float(freq[:-1]) > thresh:
                    count += 1
            if count not in counts:
                counts[count] = 0
            counts[count] += 1

        xs = []
        ys = []
        for count, num in counts.items():
            xs.append(count)
            ys.append(num)
        xs = np.array(xs)
        ys = np.array(ys)

        figure = bokeh.plotting.figure(
            title='Histogram of STR loci by number of non-rare alleles',
            x_axis_label=f'Num. alleles with frequency > {thresh*100}%',
            y_axis_label='Num. STRs analyzed',
            width=600,
            height=600,
            tools=''
        )
        figure.background_fill_color = None
        figure.border_fill_color = None
        figure.ygrid.grid_line_color = None
        figure.xgrid.grid_line_color = None
        figure.quad(
            top=ys,
            bottom=0,
            left=xs-0.5,
            right=xs+0.5
        )
        bokeh.io.export_png(figure, filename=f'{ukb}/export_scripts/results/allele_count_thresh_{thresh}.png')
        bokeh.io.export_svg(figure, filename=f'{ukb}/export_scripts/results/allele_count_thresh_{thresh}.svg')

if __name__ == '__main__':
    main()
