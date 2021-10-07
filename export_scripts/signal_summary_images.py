#!/usr/bin/env python3

import argparse
import os

import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import bokeh.transform
import matplotlib.colors
import numpy as np

ukb = os.environ['UKB']

def parse_p(f):
    if f == 'nan':
        return np.nan
    else:
        f = float(f)
        if f < 1e-100:
            f = 1e-100
        return -np.log10(f)

def get_var_peaks(peaks, var):
    return [(p_val, other_p_val)
            for peak_list in peaks.values()
            for var_type, p_val, other_p_val in peak_list
            if var_type == var and p_val != other_p_val]

def peaks_by_pheno(peaks):
    sorted_phenos = sorted(peaks, key = lambda pheno: len([peak for peak in peaks[pheno]]))
    data=dict(
        phenos=sorted_phenos,
        str_peaks=[len([peak for peak in peaks[pheno] if peak[0] == 'STR']) for pheno in sorted_phenos],
        snp_peaks=[len([peak for peak in peaks[pheno] if peak[0] == 'SNP']) for pheno in sorted_phenos]
    )

    fig = bokeh.plotting.figure(
        x_range=sorted_phenos,
        title='Peaks',
        tools='',
        x_axis_label = 'phenotype',
        y_axis_label = 'n peaks',
        width=1200,
        height=1200
    )
    fig.xaxis.major_label_orientation = 1.2
    fig.xaxis.group_label_orientation = 1.2
    fig.vbar_stack(
        [ 'snp_peaks', 'str_peaks'],
        x='phenos',
        legend_label = ['snp_peaks', 'str_peaks'],
        color=['blue', 'red'],
        source = data,
        width=0.9
    )
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.ygrid.grid_line_color = None
    fig.xgrid.grid_line_color = None

    return fig

def get_str_tagging_peak_hist(peaks):
    fig = bokeh.plotting.figure(
        title=f'Peaks',
        tools='',
        x_axis_label = '-log10 peak p-value',
        y_axis_label = 'n peaks',
        width=1200,
        height=1200
    )

    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.ygrid.grid_line_color = None
    fig.xgrid.grid_line_color = None
    peaks = sum(peaks.values(), [])

    bin_size = 2
    edges = np.arange(min(p_val for (_, p_val, _) in peaks), max(p_val for (_, p_val, _) in peaks) + bin_size, bin_size)
    bin_offset = 0.2
    peak_bins = [None]*(len(edges) - 1)

    for bin_ in range(len(peak_bins)):
        peak_bins[bin_] = []
    for var, p_val, other_p_val in peaks:
        found_bin = False
        for curr_bin, edge in enumerate(edges[1:]):
            if p_val < edge:
                bin_ = curr_bin
                found_bin = True
                break
        if not found_bin:
            bin_ = len(peak_bins) - 1
    
        if var == 'SNP':
            if np.isnan(other_p_val):
                tag = 0
            else:
                tag = other_p_val/p_val
        else:
            assert var == 'STR'
            if np.isnan(other_p_val):
                tag = np.inf
            else:
                tag = p_val/other_p_val

        peak_bins[bin_].append(tag)

    peak_bins = [np.sort(peak_bin) for peak_bin in peak_bins]
    #max_tag = min(max(tag for peak_bin in peak_bins for tag in peak_bin if ~np.isinf(tag)), 1.21)
    max_tag = 1.21
    transition = 1/max_tag
    mpl_cmap = matplotlib.colors.LinearSegmentedColormap(
        'my_cmap',
        {
         #gold
         'red':   [[0, 0, 0], [transition, 1, 1], [1, 1, 1]],
         'green': [[0, 0, 0], [transition, 0, 0], [1, 215/256, 215/256]],
         'blue':  [[0, 1, 1], [transition, 0, 0], [1, 0, 0]],
         #moccassin
         #'red':   [[0, 0, 0], [transition, 1, 1], [1, 1, 1]],
         #'green': [[0, 0, 0], [transition, 0, 0], [1, 228/256, 228/256]],
         #'blue':  [[0, 1, 1], [transition, 0, 0], [1, 181/256, 181/256]],
        }
    )
    palette=tuple('#%02x%02x%02x' % mpl_cmap(i, bytes=True)[:3] for i in range(256))

    cmap = bokeh.transform.linear_cmap('peaks', palette=palette, low=0, high=max_tag)
    first = True
    for peak_bin, edge1, edge2 in zip(peak_bins, edges[:-1], edges[1:]):
        '''
        palette= np.stack((np.arange(256), np.zeros(256), np.arange(256)[::-1]), axis=1)
        peak_colors = palette[np.floor(peak_bin*256).astype(int), :]
        peak_colors = [peak_colors[i, :] for i in range(len(peak_colors))]
        '''
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
        '''
        fig.quad(
            bottom = np.arange(0, height),
            top = np.arange(1, height+1),
            left = edge1,
            right = edge2,
            color=cmap,
            #color=peak_colors,
            legend_label=f'percent height of {alt_var} within {var} peak'
        )
        '''
        if first:
            color_bar = bokeh.models.ColorBar(
                color_mapper = cmap['transform'],
                width=8,
                title='fraction height of STR within SNP peak'
            )
            fig.add_layout(color_bar, 'right')
            first=False

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
            '''
            palette= np.stack((np.arange(256), np.zeros(256), np.arange(256)[::-1]), axis=1)
            peak_colors = palette[np.floor(peak_bin*256).astype(int), :]
            peak_colors = [peak_colors[i, :] for i in range(len(peak_colors))]
            '''
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
            '''
            fig.quad(
                bottom = np.arange(0, height),
                top = np.arange(1, height+1),
                left = edge1,
                right = edge2,
                color=cmap,
                #color=peak_colors,
                legend_label=f'percent height of {alt_var} within {var} peak'
            )
            '''
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

    '''
    sorted_peaks = sorted(peaks, key = lambda pheno: len([peak for peak in peaks[pheno] if peak[0] == 'SNP']))
    factors = [(pheno, var) for pheno in sorted_peaks for var in ('SNPs', 'STRs')]
    fig = bokeh.plotting.figure(
        x_range=bokeh.models.FactorRange(*factors),
        title='Peaks',
        tools='',
        x_axis_label = 'phenotype',
        y_axis_label = 'n peaks',
    )
    fig.xaxis.major_label_orientation = 1.2
    fig.xaxis.group_label_orientation = 1.2
    str_peaks=[len([peak for peak in peaks[pheno] if peak[0] == 'STR']) for pheno in sorted_peaks]
    snp_peaks=[len([peak for peak in peaks[pheno] if peak[0] == 'SNP']) for pheno in sorted_peaks]
    data=dict(
        peaks=[bar for pheno_bar in zip(snp_peaks, str_peaks) for bar in pheno_bar],
        factors=factors
    )
    fig.vbar(
        x='factors',
        top='peaks',
        source = data,
        width=0.9
    )
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.ygrid.grid_line_color = None
    fig.xgrid.grid_line_color = None

    figs.append(fig)
    '''

    bokeh.io.export_png(peaks_by_pheno(peaks), filename=f'{ukb}/export_scripts/results/peaks_by_pheno.png')
    bokeh.io.export_svg(peaks_by_pheno(peaks), filename=f'{ukb}/export_scripts/results/peaks_by_pheno.svg')

    str_tagging_peak_hist = get_str_tagging_peak_hist(peaks) 
    bokeh.io.export_png(str_tagging_peak_hist, filename=f'{ukb}/export_scripts/results/str_tagging_peak_hist.png')
    bokeh.io.export_svg(str_tagging_peak_hist, filename=f'{ukb}/export_scripts/results/str_tagging_peak_hist.svg')


    overall_fig = bokeh.layouts.gridplot(figs, ncols=3, plot_width = 1200, plot_height=1200, toolbar_location=None)

    bokeh.io.export_png(overall_fig, filename=f'{ukb}/export_scripts/results/matched_peak_hists.png')

if __name__ == '__main__':
    main()
