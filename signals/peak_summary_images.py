#!/usr/bin/env python3

import argparse
import math

import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import bokeh.transform
import numpy as np

import graphing_utils
import phenotypes

rng = np.random.default_rng(13)

def parse_p(f):
    if f == 'nan':
        return np.nan
    else:
        f = float(f)
        if f < 1e-300:
            f = 1e-300
        return -np.log10(f)

def get_var_peaks(peaks, var):
    return [(p_val, other_p_val)
            for peak_list in peaks.values()
            for var_type, p_val, other_p_val in peak_list
            if var_type == var and p_val != other_p_val]

def get_peaks_by_pheno(peaks):
    pheno_pairs = [('serum biomarkers', pheno) if phenotypes.is_serum_biomarker(pheno) else ('haematology', pheno) for pheno in peaks]
    sorted_pheno_pairs = sorted(pheno_pairs, key = lambda pheno_pair: (pheno_pair[0], len(peaks[pheno_pair[1]])))
    str_only_peaks=[len([peak for peak in peaks[pheno] if peak[0] == 'STR' and np.isnan(peak[2])]) for _, pheno in sorted_pheno_pairs]
    snp_only_peaks=[len([peak for peak in peaks[pheno] if peak[0] == 'SNP' and np.isnan(peak[2])]) for _, pheno in sorted_pheno_pairs]
    both_peaks=[len([peak for peak in peaks[pheno] if ~np.isnan(peak[2])]) for _, pheno in sorted_pheno_pairs]

    data=dict(
        pheno_pairs=sorted_pheno_pairs,
        str_only_peaks=str_only_peaks,
        snp_only_peaks=snp_only_peaks,
        both_peaks=both_peaks
    )

    fig = bokeh.plotting.figure(
        x_range=bokeh.models.FactorRange(*sorted_pheno_pairs),
        title='Peaks',
        tools='',
        x_axis_label = 'phenotype',
        y_axis_label = 'n peaks',
        width=math.floor(4.25/2*1200),
        height=1200,
    )
    fig.toolbar_location = None
    fig.xaxis.major_label_orientation = 1.2
    fig.x_range.range_padding = 0.1
    fig.vbar_stack(
        [ 'snp_only_peaks', 'both_peaks', 'str_peaks'],
        x='pheno_pairs',
        legend_label = ['Only tagged by SNPs', 'Tagged by both SNPs and STRs', 'Only tagged by STRs'],
        color=['#00B8FF', '#9D6898', '#FF520D'],
        source = data,
        width=0.9
    )
    fig.legend.location='top_left'
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color = None
    fig.title.text_font_size = '30px'
    fig.axis.axis_label_text_font_size = '26px'
    fig.axis.major_label_text_font_size = '20px'
    fig.xaxis.group_text_font_size = '26px'
    fig.legend.label_text_font_size = '18px'

    graphing_utils.resize(fig, 5000/1200)
    return fig

def get_p_val_comp_scatter(peaks):
    fig = bokeh.plotting.figure(
        title='STR vs SNP/indel peak p-values',
        x_axis_label = '-log10 SNP/indel peak p-values',
        y_axis_label = '-log10 STR peak p-values',
        width=1200,
        height=1200
    )

    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color = None
    fig.toolbar_location = None
    fig.title.text_font_size = '30px'
    fig.axis.axis_label_text_font_size = '26px'
    fig.axis.major_label_text_font_size = '20px'
    peaks = sum(peaks.values(), [])

    snp_ps = []
    str_ps = []
    for var_type, p1, p2 in peaks:
        if np.isnan(p1):
            p1 = 0
        if np.isnan(p2):
            p2 = 0
        if var_type == 'SNP':
            snp_ps.append(p1)
            str_ps.append(p2)
        else:
            assert var_type == 'STR'
            snp_ps.append(p2)
            str_ps.append(p1)

    snp_ps = np.array(snp_ps)
    str_ps = np.array(str_ps)

    maxed = (snp_ps == 300) & (str_ps == 300)
    snp_ps[maxed] += rng.standard_normal(np.sum(maxed))*9
    str_ps[maxed] += rng.standard_normal(np.sum(maxed))*9

    fig.circle(snp_ps, str_ps, radius=3, alpha=.5, color='#d2691e')
    graphing_utils.resize(fig, 5000/1200, legend=False)
    return fig

def linear_int_interpolate(c1, c2, dist):
    c_new = []
    for coord1, coord2 in zip(c1, c2):
        c_new.append(coord1 + round((coord2 - coord1)*dist))
    return c_new

def get_p_val_heatmap(peaks):
    fig = bokeh.plotting.figure(
        title='STR vs SNP/indel peak p-values',
        x_axis_label = '-log10 SNP/indel peak p-values',
        y_axis_label = '-log10 STR peak p-values',
        width=1200,
        height=1200,
        output_backend='svg'
    )

    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color = None
    fig.toolbar_location = None
    fig.title.text_font_size = '30px'
    fig.axis.axis_label_text_font_size = '26px'
    fig.axis.major_label_text_font_size = '20px'
    peaks = sum(peaks.values(), [])

    snp_ps = []
    str_ps = []

    for var_type, p1, p2 in peaks:
        if np.isnan(p1):
            p1 = 0
        if np.isnan(p2):
            p2 = 0
        if var_type == 'SNP':
            snp_ps.append(p1)
            str_ps.append(p2)
        else:
            assert var_type == 'STR'
            snp_ps.append(p2)
            str_ps.append(p1)

    snp_ps = np.array(snp_ps)
    str_ps = np.array(str_ps)

    bin_size = 8
    max_p = 296
    grid = np.mgrid[:(max_p+bin_size):bin_size, :(max_p+bin_size):bin_size]
    x = grid[0,:,:].flatten()
    y = grid[1,:,:].flatten()
    x, y = x[(x >= bin_size) | (y >= bin_size)], y[(x >= bin_size) | (y >= bin_size)]

    cds = bokeh.models.ColumnDataSource(dict(
        left = x,
        right = x + bin_size,
        bottom = y,
        top = y + bin_size,
        counts = np.sum(
            (snp_ps[:, np.newaxis] >= x[np.newaxis, :]) &
            (snp_ps[:, np.newaxis] < x[np.newaxis, :] + bin_size) &
            (str_ps[:, np.newaxis] >= y[np.newaxis, :]) &
            (str_ps[:, np.newaxis] < y[np.newaxis, :] + bin_size),
        axis = 0)
    ))
    palette = [
        linear_int_interpolate((134,204,195), (9,41,46), i/254) for i in range(-1, 255)
    ]
    cmap = bokeh.transform.log_cmap(
        'counts',
        palette = palette,
        low=1,
        high=np.max(cds.data['counts']),
        low_color=(255, 255, 255)
    )
    color_mapper = bokeh.models.LogColorMapper(
        palette = palette,
        low=1,
        high=np.max(cds.data['counts'])
    )

    fig.quad(
        left='left', right='right', bottom='bottom', top='top', source=cds, fill_color=cmap, line_width=0
    )
    fig.line(
        x = [0, max_p + bin_size],
        y = [0, max_p + bin_size],
        line_width = 20,
        color = 'dimgrey',
    )

    color_bar = bokeh.models.ColorBar(
        color_mapper = color_mapper,
        width=70,
        major_label_text_font_size = '20px'
    )
    fig.add_layout(color_bar, 'right')

    graphing_utils.resize(fig, 5000/1200, legend=False)
    return fig

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('--phenotypes', nargs='+')
    parser.add_argument('--peak-files', nargs='+')
    args = parser.parse_args()
    assert len(args.phenotypes) == len(args.peak_files)

    peaks = {}
    for phenotype, peak_fname in zip(args.phenotypes, args.peak_files):
        peaks[phenotype] = []
        with open(peak_fname) as peaks_file:
            header = next(peaks_file).strip().split('\t')
            var_type_idx = header.index('variant_type')
            p_val_idx = header.index('p_value')
            p_val_other_idx = header.index('p_value_other_variant_type')

            for line in peaks_file:
                split = line.strip().split('\t')
                peaks[phenotype].append((
                    split[var_type_idx],
                    parse_p(split[p_val_idx]),
                    parse_p(split[p_val_other_idx])
                ))

    peaks_by_pheno = get_peaks_by_pheno(peaks)
    bokeh.io.export_png(peaks_by_pheno, filename=f'{args.outdir}/peaks_by_pheno.png')
    bokeh.io.export_svg(peaks_by_pheno, filename=f'{args.outdir}/peaks_by_pheno.svg')

    p_val_heatmap = get_p_val_heatmap(peaks)
    bokeh.io.export_png(p_val_heatmap, filename=f'{args.outdir}/peak_p_val_heatmap.png')
    bokeh.io.export_svg(p_val_heatmap, filename=f'{args.outdir}/peak_p_val_heatmap.svg')

if __name__ == '__main__':
    main()
