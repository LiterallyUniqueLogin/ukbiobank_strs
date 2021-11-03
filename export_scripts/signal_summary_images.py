#!/usr/bin/env python3

import argparse
import math
import os

import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import bokeh.transform
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

import graphing_utils

ukb = os.environ['UKB']
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
    sorted_phenos = sorted(peaks, key = lambda pheno: len(peaks[pheno]))
    str_only_peaks=[len([peak for peak in peaks[pheno] if peak[0] == 'STR' and np.isnan(peak[2])]) for pheno in sorted_phenos]
    snp_only_peaks=[len([peak for peak in peaks[pheno] if peak[0] == 'SNP' and np.isnan(peak[2])]) for pheno in sorted_phenos]
    both_peaks=[len([peak for peak in peaks[pheno] if ~np.isnan(peak[2])]) for pheno in sorted_phenos]

    percent_str_only = sum(str_only_peaks)/(sum(snp_only_peaks) + sum(both_peaks) + sum(str_only_peaks))
    percent_snp_only = sum(snp_only_peaks)/(sum(snp_only_peaks) + sum(both_peaks) + sum(str_only_peaks))

    data=dict(
        phenos=sorted_phenos,
        str_only_peaks=str_only_peaks,
        snp_only_peaks=snp_only_peaks,
        both_peaks=both_peaks
    )

    fig = bokeh.plotting.figure(
        x_range=sorted_phenos,
        title='Peaks',
        tools='',
        x_axis_label = 'phenotype',
        y_axis_label = 'n peaks',
        width=math.floor(4.25/2*1200),
        height=1200,
    )
    fig.toolbar_location = None
    fig.xaxis.major_label_orientation = 1.2
    fig.xaxis.group_label_orientation = 1.2
    fig.vbar_stack(
        [ 'snp_only_peaks', 'both_peaks', 'str_peaks'],
        x='phenos',
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotypes', nargs='+')
    phenotypes = parser.parse_args().phenotypes

    peaks = {}
    for phenotype in phenotypes:
        peaks[phenotype] = []
        with open(f'{ukb}/signals/peaks/{phenotype}_250000_5e-8.tab') as peaks_file:
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

    figs = []

    for var, alt_var in (('SNP', 'STR'), ('STR', 'SNP')):
        fig = bokeh.plotting.figure(
            title=f'{var} peaks',
            tools='',
            x_axis_label = '-log10 peak significance',
            y_axis_label = 'n peaks'
        )

        fig.background_fill_color = None
        fig.border_fill_color = None
        fig.ygrid.grid_line_color = None
        fig.xgrid.grid_line_color = None

        var_peaks = get_var_peaks(peaks, var)
        matched_var_peaks = [p_val for p_val, other_p_val in var_peaks if ~np.isnan(other_p_val)]
        unmatched_var_peaks = [p_val for p_val, other_p_val in var_peaks if np.isnan(other_p_val)]

        bin_size = 2
        edges = np.arange(min(p_val for (p_val, _) in var_peaks), max(p_val for (p_val, _) in var_peaks) + bin_size, bin_size)
        bin_offset = 0.2

        matched_hist, _ = np.histogram(matched_var_peaks, edges)
        unmatched_hist, _ = np.histogram(unmatched_var_peaks, edges)
        centers = (edges[:-1] + edges[1:])/2

        fig.quad(
            left=edges[:-1]+bin_offset,
            bottom=matched_hist,
            right=edges[1:]-bin_offset,
            top=unmatched_hist+matched_hist,
            color='blue',
            legend_label=f'no sig {alt_var} in peak'
        )
        fig.quad(
            left=edges[:-1]+bin_offset,
            bottom=np.zeros((len(edges)-1, )),
            right=edges[1:]-bin_offset,
            top=matched_hist,
            color='red',
            legend_label=f'sig {alt_var} in peak'
        )
        figs.append(fig)

        fig = bokeh.plotting.figure(
            title=f'{var} peaks',
            tools='',
            x_axis_label = '-log10 peak significance',
            y_axis_label = 'n peaks'
        )

        fig.background_fill_color = None
        fig.border_fill_color = None
        fig.ygrid.grid_line_color = None
        fig.xgrid.grid_line_color = None

        fig.quad(
            top=unmatched_hist + matched_hist,
            left=edges[:-1]+bin_offset,
            right=edges[1:]-bin_offset,
            bottom=np.zeros((len(edges) - 1, )),
            color='blue',
            legend_label = f'all {var} peaks'
        )

        peak_bins = [None]*(len(edges) - 1)
        for bin_ in range(len(peak_bins)):
            peak_bins[bin_] = []
        for p_val, other_p_val in var_peaks:
            found_bin = False
            for curr_bin, edge in enumerate(edges[1:]):
                if p_val < edge:
                    bin_ = curr_bin
                    found_bin = True
                    break
            if not found_bin:
                bin_ = len(peak_bins) - 1
            if ~np.isnan(other_p_val):
                peak_bins[bin_].append(other_p_val/p_val)
            else:
                peak_bins[bin_].append(0)

        peak_bins = [np.sort(peak_bin)[::-1] for peak_bin in peak_bins]

        for peak_bin, edge1, edge2 in zip(peak_bins, edges[:-1], edges[1:]):
            height = len(peak_bin)
            n_zeros = sum(frac == 0 for frac in peak_bin)
            if height-n_zeros > 0:
                edge2 -= bin_offset
                edge1 += bin_offset
                diff = (edge2 - edge1)/2
                fig.quad(
                    bottom = np.arange(0, height-n_zeros),
                    top = np.arange(1, height-n_zeros+1),
                    left = edge1+diff*(1-peak_bin[:(height-n_zeros)]),
                    right = edge2-diff*(1-peak_bin[:(height-n_zeros)]),
                    color='red',
                    legend_label=f'percent height of {alt_var} within {var} peak'
                )

        figs.append(fig)

        fig = bokeh.plotting.figure(
            title=f'{var} peaks',
            tools='',
            x_axis_label = '-log10 peak significance',
            y_axis_label = 'n peaks'
        )

        fig.background_fill_color = None
        fig.border_fill_color = None
        fig.ygrid.grid_line_color = None
        fig.xgrid.grid_line_color = None

        peak_bins = [None]*(len(edges) - 1)
        for bin_ in range(len(peak_bins)):
            peak_bins[bin_] = []
        for p_val, other_p_val in var_peaks:
            found_bin = False
            for curr_bin, edge in enumerate(edges[1:]):
                if p_val < edge:
                    bin_ = curr_bin
                    found_bin = True
                    break
            if not found_bin:
                bin_ = len(peak_bins) - 1
            if ~np.isnan(other_p_val):
                peak_bins[bin_].append(other_p_val/p_val)
            else:
                peak_bins[bin_].append(0)

        peak_bins = [np.sort(peak_bin)[::-1] for peak_bin in peak_bins]
        cmap = bokeh.transform.linear_cmap('peaks', palette=tuple('#%02x%02x%02x' % (i, 0, 255-i) for i in range(256)), low=0, high=1)
        first = True
        for peak_bin, edge1, edge2 in zip(peak_bins, edges[:-1], edges[1:]):
            height = len(peak_bin)
            edge2 -= bin_offset
            edge1 += bin_offset
            data=dict(
                bottom = np.arange(0, height),
                top = np.arange(1, height+1),
                left = np.full(height, edge1),
                right = np.full(height, edge2),
                peaks=peak_bin
            )
            fig.quad(
                bottom = 'bottom',
                top = 'top',
                left = 'left',
                right = 'right',
                color=cmap,
                source=data
            )

            if first:
                color_bar = bokeh.models.ColorBar(
                    color_mapper = cmap['transform'],
                    width=8,
                    title='fraction height of STR within SNP peak'
                )
                fig.add_layout(color_bar, 'right')
                first=False

        figs.append(fig)

    sorted_peaks = sorted(peaks, key = lambda pheno: len([peak for peak in peaks[pheno] if peak[0] == 'SNP']))
    fig = bokeh.plotting.figure(
        x_range=sorted_peaks,
        title='SNP peaks',
        tools='',
        x_axis_label = 'phenotype',
        y_axis_label = 'n peaks',
    )
    fig.xaxis.major_label_orientation = 1.2
    data=dict(
        phenos=sorted_peaks,
        unmatched_peaks=[len([peak for peak in peaks[pheno] if peak[0] == 'SNP' and np.isnan(peak[2])]) for pheno in sorted_peaks],
        matched_peaks=[len([peak for peak in peaks[pheno] if peak[0] == 'SNP' and ~np.isnan(peak[2])]) for pheno in sorted_peaks],
    )
    fig.vbar_stack(
        ['matched_peaks','unmatched_peaks'],
        x='phenos',
        legend_label = ['matched_peaks', 'unmatched_peaks'],
        color=['red', 'blue'],
        source = data,
        width=0.9
    )
    
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.ygrid.grid_line_color = None
    fig.xgrid.grid_line_color = None

    figs.append(fig)

    peaks_by_pheno = get_peaks_by_pheno(peaks)
    bokeh.io.export_png(peaks_by_pheno, filename=f'{ukb}/export_scripts/results/peaks_by_pheno.png')
    bokeh.io.export_svg(peaks_by_pheno, filename=f'{ukb}/export_scripts/results/peaks_by_pheno.svg')

    p_val_comp_scatter = get_p_val_comp_scatter(peaks)
    bokeh.io.export_png(p_val_comp_scatter, filename=f'{ukb}/export_scripts/results/p_val_comp_jitter_scatter.png')
    bokeh.io.export_svg(p_val_comp_scatter, filename=f'{ukb}/export_scripts/results/p_val_comp_jitter_scatter.svg')

if __name__ == '__main__':
    main()
