#!/usr/bin/env python3

import argparse
import math
import time

import bokeh.io
import bokeh.models
import bokeh.plotting
import numpy as np
import pandas as pd
import polars as pl

import graphing_utils
import region_plots

overview_max_p_val = 100

def make_overview_manhattan(
    outfname,
    phenotype,
    my_str_results,
    plink_snp_results,
    legendless,
    chr_lens,
    has_peaks,
    maf_threshold,
    snp_maf_file,
):
    print(f"Plotting overview of phenotype {phenotype} ... ", end='', flush=True)
    start_time = time.time()

    my_str_results['p_val'] = np.minimum(overview_max_p_val, my_str_results['p_val'])
    plink_snp_results['p_val'] = np.minimum(overview_max_p_val, plink_snp_results['p_val'])

    if has_peaks:
        str_peak_dict = {}
        snp_peak_dict = {}
    str_locus_dict = {}
    snp_locus_dict = {}
    insignificant_locus_dict = {}

    copy_cols = {'chr', 'pos', 'p_val'}
    gwas_sig = -np.log10(5e-8)

    for col in copy_cols:
        if has_peaks:
            str_peak_dict[col] = my_str_results[my_str_results['is_peak']][col]
            snp_peak_dict[col] = plink_snp_results[plink_snp_results['is_peak']][col]

        str_locus_dict[col] = my_str_results[my_str_results['p_val'] > gwas_sig][col]
        snp_locus_dict[col] = plink_snp_results[plink_snp_results['p_val'] > gwas_sig][col]

        insignificant_locus_dict[col] = np.concatenate((
            my_str_results[my_str_results['p_val'] < gwas_sig][col],
            plink_snp_results[plink_snp_results['p_val'] < gwas_sig][col]
        ))

    for _dict, color1, color2 in [
            (insignificant_locus_dict, '#A9A9A9', '#696969')]:
        _dict['color'] = np.full(_dict['chr'].shape, '#000000')
        _dict['color'][_dict['chr'] % 2 == 1] = color1
        _dict['color'][_dict['chr'] % 2 == 0] = color2

    manhattan_plot = bokeh.plotting.figure(
        width=math.floor(4.25*400),
        height=400,
        title=(phenotype.capitalize()),
        x_axis_label='Chromosomes',
        y_axis_label='-log10(p-value)',
        output_backend='canvas' if outfname[-4:] == '.png' else 'svg'
    )
    manhattan_plot.title.text_font_size = '30px'
    manhattan_plot.axis.axis_label_text_font_size = '26px'
    manhattan_plot.axis.major_label_text_font_size = '20px'
    manhattan_plot.grid.grid_line_color = None
    manhattan_plot.background_fill_color = None
    manhattan_plot.border_fill_color = None
    manhattan_plot.toolbar_location = None

    region_plots.full_genome_x_axis(
        manhattan_plot, chr_lens
    )
    dict_list = [str_locus_dict, snp_locus_dict, insignificant_locus_dict]
    if has_peaks:
        dict_list += [str_peak_dict, snp_peak_dict]
    for dict_ in dict_list:
        region_plots.full_genome_pandas_df(dict_, chr_lens)

    # from https://projects.susielu.com/viz-palette
    str_color = '#FF520D'
    snp_color = '#00B8FF'
    size_ratio = 5000/1200
    
    manhattan_plot.diamond(
        'plot_pos',
        'p_val',
        source=bokeh.models.ColumnDataSource(snp_locus_dict),
        color=snp_color,
        size=math.floor(size_ratio*8)
    )
    manhattan_plot.diamond(
        'plot_pos',
        'p_val',
        source=bokeh.models.ColumnDataSource(str_locus_dict),
        color=str_color,
        size=math.floor(size_ratio*8),
        alpha=0.5
    )
    if has_peaks:
        if not legendless:
            manhattan_plot.diamond(
                'plot_pos',
                'p_val',
                source=bokeh.models.ColumnDataSource(snp_peak_dict),
                fill_color=snp_color,
                legend_label = 'SNPs',
                size=math.floor(size_ratio*15)
            )
        else:
            manhattan_plot.diamond(
                'plot_pos',
                'p_val',
                source=bokeh.models.ColumnDataSource(snp_peak_dict),
                fill_color=snp_color,
                size=math.floor(size_ratio*15)
            )
        if not legendless:
            manhattan_plot.diamond(
                'plot_pos',
                'p_val',
                source=bokeh.models.ColumnDataSource(str_peak_dict),
                fill_color=str_color,
                legend_label = 'STRs',
                size=math.floor(size_ratio*15)
            )
        else:
            manhattan_plot.diamond(
                'plot_pos',
                'p_val',
                source=bokeh.models.ColumnDataSource(str_peak_dict),
                fill_color=str_color,
                size=math.floor(size_ratio*15)
            )
    if not legendless:
        manhattan_plot.diamond(
            'plot_pos',
            'p_val',
            source=bokeh.models.ColumnDataSource(insignificant_locus_dict),
            color='color',
            size=math.floor(size_ratio*8),
            legend_label = 'not GWAS sig.'
        )
    else:
        manhattan_plot.diamond(
            'plot_pos',
            'p_val',
            source=bokeh.models.ColumnDataSource(insignificant_locus_dict),
            color='color',
            size=math.floor(size_ratio*8),
        )

    if not legendless:
        manhattan_plot.legend.label_text_font_size = '22px'

    graphing_utils.resize(manhattan_plot, size_ratio, legend=not legendless)

    if outfname[-4:] == '.png':
        bokeh.io.export_png(manhattan_plot, filename=outfname)
    else:
        assert outfname[-4:] == '.svg'
        bokeh.io.export_svg(manhattan_plot, filename=outfname)

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('outfile')
    parser.add_argument('my_str_results')
    parser.add_argument('plink_snp_results')
    parser.add_argument('peaks_fname')
    parser.add_argument('chr_lens')
    parser.add_argument('--legendless', action='store_true', default=False)
    parser.add_argument('--binary', default=False, choices={'logistic', 'firth'})
    parser.add_argument('--maf-threshold', type=float)
    parser.add_argument('--snp-maf-file')
    parser.add_argument('--info-threshold', type=float)
    parser.add_argument('--snp-info-files', nargs=22)
    args = parser.parse_args()

    peaks_fname = args.peaks_fname
    if peaks_fname == '.':
        peaks_fname = None

    assert (args.snp_maf_file is None) == (args.maf_threshold is None)

    phenotype = args.phenotype
    binary = args.binary

    my_str_results = region_plots.load_my_str_results(
        phenotype, binary, args.my_str_results
    )
    plink_snp_results = region_plots.load_plink_results(
        phenotype, binary, args.plink_snp_results
    )

    if peaks_fname is not None:
        my_str_results, plink_snp_results = region_plots.load_and_merge_peaks_into_dfs(
            peaks_fname,
            my_str_results,
            plink_snp_results
        )

    if not isinstance(plink_snp_results, pl.DataFrame):
        plink_snp_results = pl.DataFrame(plink_snp_results)

    if args.maf_threshold is not None:
        mafs = pl.scan_csv(
            args.snp_maf_file,
            separator='\t',
        )
        plink_snp_results = plink_snp_results.lazy().join(
            mafs,
            how='left',
            left_on=['chr', 'pos', 'ref', 'alt'],
            right_on=['#CHROM', 'POS', 'REF', 'ALT']
        ).filter(
            pl.col('ALT_FREQS') > args.maf_threshold
        ).collect()

    if args.info_threshold is not None:
        infos_per_chrom = [pl.read_csv(
            info_file,
            has_header = False,
            new_columns = ['combined_ID', 'ID', 'pos', 'ref', 'alt', 'maf', '?', 'info'],
            separator='\t',
            null_values='NA'
        ).with_columns(pl.lit(chrom).alias('chr').cast(int)) for (chrom, info_file) in zip(range(1, 23), args.snp_info_files)]

        infos = infos_per_chrom[0]
        for other_infos in infos_per_chrom[1:]:
            infos = infos.vstack(other_infos)

        plink_snp_results = plink_snp_results.join(
			infos,
            how='left',
            on=['chr', 'pos', 'ref', 'alt'],
        ).filter(
            pl.col('info') > args.info_threshold
        )

    plink_snp_results = plink_snp_results.to_numpy(structured=True)

    chr_lens = np.genfromtxt(
        args.chr_lens,
        usecols=[1],
        skip_header=1,
        dtype=int
    )

    make_overview_manhattan(
        args.outfile,
        phenotype,
        my_str_results,
        plink_snp_results,
        args.legendless,
        chr_lens,
        peaks_fname is not None,
        args.maf_threshold,
        args.snp_maf_file,
    )

if __name__ == "__main__":
    main()

