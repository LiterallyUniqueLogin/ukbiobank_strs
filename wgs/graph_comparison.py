#!/usr/bin/env python3

import argparse

import bokeh.io
import bokeh.plotting
import numpy as np
import polars as pl

def compare_calls(args):
    df = pl.read_csv(args.locus_compare, sep='\t')
    n_loci = df.shape[0]

    finemapped_df = pl.read_csv(args.conf_finemapped_strs, sep='\t').select([
        'chrom',
        'start_pos',
        pl.lit(True).alias('finemapped')
    ]).unique()

    df = df.join(
        finemapped_df,
        how='left',
        left_on=['chrom', 'start'],
        right_on=['chrom', 'start_pos']
    )

    for stat_name in [
        'fraction_concordant_len',
        'fraction_concordant_len_sum',
        'balanced_accuracy',
        'mean_absolute_difference',
        'r',
        'dosage_r',
        'numcalls'
    ]:
        fig = bokeh.plotting.figure(
            title=f'{stat_name} histogram',
            y_axis_label = f'number of loci (out of {n_loci})',
            x_axis_label = stat_name,
            tools='',
            toolbar_location = None,
            width=800,
            height=500,
            output_backend='svg'
        )

        if stat_name != 'numcalls':
            bin_size = 0.01
        else:
            bin_size = 5000

        stat = df[stat_name].to_numpy()
        min_bin_size_multiplier = np.min(stat // bin_size)
        max_bin_size_multiplier = np.max(stat // bin_size)
        bins = np.arange(min_bin_size_multiplier, max_bin_size_multiplier + 2) * bin_size
        prev_ys = np.zeros(len(bins) - 1)

        for (finemapped, color, legend) in (True, 'teal', 'confidently fine-mapped'), (False, 'blue', 'singly fine-mapped'):
            curr_stat = df.filter(pl.col('finemapped').is_null() != finemapped)[stat_name].to_numpy()
            ys = np.sum((curr_stat.reshape(-1, 1) >= bins[:-1]) & (curr_stat.reshape(-1, 1) < bins[1:]), axis=0)

            fig.quad(top=(ys + prev_ys)[ys != 0], bottom=prev_ys[ys != 0], left=bins[:-1][ys != 0], right=bins[1:][ys != 0], color=color, legend_label=legend, line_width=.9)
            prev_ys += ys

        bokeh.io.export_png(fig, filename=f'{args.out}/{stat_name}.png')
        bokeh.io.export_svg(fig, filename=f'{args.out}/{stat_name}.svg')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--locus-compare', default='wgs/compare_calls_w_stage-0/out-locuscompare.tab')
    parser.add_argument('--conf-finemapped-strs', default='post_finemapping/results/confidently_finemapped_strs_for_paper.tab')
    parser.add_argument('--out', default='wgs/compare_calls_w_stage-0')
    args = parser.parse_args()

    compare_calls(args)

if __name__ == '__main__':
    main()
