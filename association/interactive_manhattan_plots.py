#!/usr/bin/env python3

import argparse
import ast
import copy
import os
import os.path
import time
from typing import Dict, Tuple, Optional, Set

import bokeh.colors
import bokeh.colors.named
import bokeh.embed
import bokeh.events
import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.models.callbacks
import bokeh.models.tools
import bokeh.plotting
import bokeh.resources
import colorcet
import numpy as np
import numpy.lib.recfunctions
import numpy.ma
import numpy.random
import pandas as pd

import python_array_utils as utils

ukb = os.environ['UKB']

chr_lens = np.genfromtxt(
    f'{ukb}/misc_data/genome/chr_lens.txt',
    skip_header=1,
    usecols=(1),
    dtype=int
)

overview_max_p_val = 100

def make_overview_manhattan(
        outfname,
        phenotype,
        my_str_results,
        plink_snp_results):
    print(f"Plotting overview of phenotype {phenotype} ... ", end='', flush=True)
    start_time = time.time()

    my_str_results['p_val'] = np.minimum(overview_max_p_val, my_str_results['p_val'])
    plink_snp_results['p_val'] = np.minimum(overview_max_p_val, plink_snp_results['p_val'])
  
    '''
    n_str_peaks = np.sum(my_str_results['is_peak'])
    n_snp_peaks = np.sum(plink_snp_results['is_peak'])
    n_strs = np.sum(~my_str_results['is_peak'])
    n_snps = np.sum(~plink_snp_results['is_peak'])
    '''

    str_peak_dict = {}
    snp_peak_dict = {}
    str_locus_dict = {}
    snp_locus_dict = {}
    insignificant_locus_dict = {}

    copy_cols = {'chr', 'pos', 'p_val'}
    gwas_sig = -np.log10(5e-8)

    for col in copy_cols:
        str_peak_dict[col] = my_str_results[my_str_results['is_peak']][col]
        snp_peak_dict[col] = plink_snp_results[plink_snp_results['is_peak']][col]

        str_locus_dict[col] = my_str_results[my_str_results['p_val'] > gwas_sig][col]
        snp_locus_dict[col] = plink_snp_results[plink_snp_results['p_val'] > gwas_sig][col]

        insignificant_locus_dict[col] = np.concatenate((
            my_str_results[my_str_results['p_val'] < gwas_sig][col],
            plink_snp_results[plink_snp_results['p_val'] < gwas_sig][col]
        ))

    for _dict, color1, color2 in [
            #(str_peak_dict,  '#EE82EE', '#DA70D6'),
            #(str_locus_dict, '#EE82EE', '#DA70D6'),
            #(snp_peak_dict,  '#87CEFA', '#1E90FF', '#6A5ACD'),
            #(snp_locus_dict, '#87CEFA', '#1E90FF', '#6A5ACD'),
            (insignificant_locus_dict, '#A9A9A9', '#696969')]:
        _dict['color'] = np.full(_dict['chr'].shape, '#000000')
        _dict['color'][_dict['chr'] % 2 == 1] = color1
        _dict['color'][_dict['chr'] % 2 == 0] = color2

    '''
    my_str_results = my_str_results[my_str_results['p_val'] > -np.log10(5e-8)]
    plink_snp_results = plink_snp_results[plink_snp_results['p_val'] > -np.log10(5e-8)]
    '''

    for chrom in range(2, 23):
        for _dict in str_peak_dict, snp_peak_dict, str_locus_dict, snp_locus_dict, insignificant_locus_dict:
            _dict['pos'][_dict['chr'] >= chrom] += chr_lens[chrom - 2]

    '''
    col_dict['type'] = np.full((n_strs + n_snps, ), 'a very long string')
    col_dict['type'][n_snps:][my_str_results['is_peak']] = 'STR peak'
    col_dict['type'][n_snps:][~my_str_results['is_peak']] = 'STR'
    col_dict['type'][:n_snps][plink_snp_results['is_peak']] = 'SNP peak'
    col_dict['type'][:n_snps][~plink_snp_results['is_peak']] = 'SNP'
    #col_dict['type'][col_dict['p_val'] < -np.log(5e-8)]  = 'Insignificant'

    col_dict['size'] = np.full((n_strs + n_snps, ), 4)
    col_dict['size'][n_snps:][my_str_results['is_peak']] = 15
    col_dict['size'][:n_snps][plink_snp_results['is_peak']] = 15

    col_dict['alpha'] = np.full((n_strs + n_snps, ), np.nan)
    col_dict['alpha'][:n_snps] = 1
    col_dict['alpha'][n_snps:] = 0.3
    col_dict['alpha'][n_snps:][my_str_results['is_peak']] = 1
    col_dict['alpha'][:n_snps][plink_snp_results['is_peak']] = 1
    '''

    manhattan_plot = bokeh.plotting.figure(
        width=1200,
        height=400,
        title=(phenotype.capitalize()),
        x_axis_label='Chromosomes',
        y_axis_label='-log10(p-value)',
        #y_range=(np.min(insignificant_locus_dict['p_val']), max(np.max(str_locus_dict['p_val']), np.max(snp_locus_dict['p_val']))),
    )
    manhattan_plot.xgrid.grid_line_color = None
    manhattan_plot.ygrid.grid_line_color = None
    manhattan_plot.background_fill_color = None
    manhattan_plot.border_fill_color = None

    pre_chr_sums = np.cumsum([0, *chr_lens[:-1]])
    mid_points = [int(num) for num in pre_chr_sums + (chr_lens//2)]
    manhattan_plot.xaxis.ticker = mid_points
    manhattan_plot.xaxis.major_label_overrides = {
        mid_points[chrom - 1]: str(chrom) for chrom in range(1, 23)
        #mid_points[chrom - 1]: f"chr{chrom}" for chrom in range(1, 23)
    }
    #manhattan_plot.xaxis.major_label_orientation = np.pi/4

    '''
    color_dict = {
        'STR peak': bokeh.colors.RGB(204, 121, 167),
        'SNP peak': bokeh.colors.RGB(86, 180, 233),
        'STR': bokeh.colors.RGB(204, 121, 167),
        'SNP': bokeh.colors.RGB(86, 180, 233),
        #'Insignificant': bokeh.colors.named.darkgrey
    }

    # from the periodic table example
    # https://docs.bokeh.org/en/latest/docs/user_guide/categorical.html#mixed-factors
    cmap = bokeh.transform.factor_cmap(
        'type',
        palette = list(color_dict.values()),
        factors = list(color_dict.keys())
    )
    '''

    # from https://projects.susielu.com/viz-palette
    str_color = "#FFB14E"
    #snp_color = "#EA5F94"
    snp_color = "#0000FF"
    manhattan_plot.diamond(
        'pos',
        'p_val',
        source=bokeh.models.ColumnDataSource(insignificant_locus_dict),
        color='color',
        size=4
    )
    manhattan_plot.diamond(
        'pos',
        'p_val',
        source=bokeh.models.ColumnDataSource(snp_locus_dict),
        color=snp_color,
        size=4
    )
    manhattan_plot.diamond(
        'pos',
        'p_val',
        source=bokeh.models.ColumnDataSource(str_locus_dict),
        color=str_color,
        size=4,
        alpha=0.5,
    )
    manhattan_plot.diamond(
        'pos',
        'p_val',
        source=bokeh.models.ColumnDataSource(snp_peak_dict),
        fill_color=snp_color,
        #line_width=2,
        legend_label = 'SNPs',
        size=11
    )
    manhattan_plot.diamond(
        'pos',
        'p_val',
        source=bokeh.models.ColumnDataSource(str_peak_dict),
        fill_color=str_color,
        #line_width=2,
        legend_label = 'STRs',
        size=11
    )

    bokeh.io.export_png(manhattan_plot, filename=outfname)

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)


# python stops being able to distinguish between small numbers and zero at 1e-324 == 0
# so cutoff a bit before then
max_p_val = 300 # in -log10

def get_conditioned_strs(condition):
    splits = condition.split('_')
    return [int(STR) for STR in splits[(splits.index('STR')+1):(splits.index('ISNP')-1)]]

def get_conditioned_isnps(condition):
    splits = condition.split('_')
    return splits[(splits.index('ISNP')+1):(splits.index('ASNP')-1)]

def create_source_dict(
        data: np.recarray,
        start_chrom: Optional[int],
        cols_to_skip: Set[str] = set(),
        cols_to_include: Set[str] = set(),
        chrs_to_var_signals = None,
        start_p_val_cap = None
    ) -> Tuple[bokeh.models.ColumnDataSource, Dict[int, bokeh.models.ColumnDataSource]]:
    if not start_chrom:
        cds_range = range(1, 23)
        start_chrom = 1
    else:
        cds_range = [start_chrom]

    sources = {}
    assert len(cols_to_skip) == 0 or len(cols_to_include) == 0
    for field in 'chr', 'pos', 'p_val':
        if field not in data.dtype.names:
            print(field, flush=True)
            assert False
    for chrom in cds_range:
        chrom_dict = {}
        idx = data['chr'] == chrom
        for name in data.dtype.names:
            if len(cols_to_skip) > 0:
                if name in cols_to_skip:
                    continue
            elif len(cols_to_include) > 0:
                if name not in cols_to_include:
                    continue
            chrom_dict[name] = data[name][idx]
        sources[chrom] = bokeh.models.ColumnDataSource(chrom_dict)
        if start_p_val_cap is not None:
            sources[chrom].data['display_p_val'] = np.minimum(
                sources[chrom].data['p_val'], start_p_val_cap
            )
        else:
            sources[chrom].data['display_p_val'] = sources[chrom].data['p_val'].copy()

        if chrs_to_var_signals:
            colnames_to_merge_on = list(chrs_to_var_signals[chrom].columns)[:-1]
            merge_cols = pd.DataFrame(data[colnames_to_merge_on][idx])

            merged_data = merge_cols.merge(
                chrs_to_var_signals[chrom].astype({'pos': int}),
                on=colnames_to_merge_on,
                how='left'
            )
            assert merged_data.shape[1] == len(chrs_to_var_signals[chrom].columns)
            assert merged_data.shape[0] == merge_cols.shape[0]
            sources[chrom].data['FINEMAP_pcausal'] = merged_data['pcausal']

    copy_source = bokeh.models.ColumnDataSource(copy.deepcopy(sources[start_chrom].data))

    return copy_source, sources

def make_manhattan_plots(
        outfname,
        phenotype,
        binary,
        unit,
        my_str_data,
        my_str_run_date,
        plink_snp_data,
        plink_snp_run_date,
        gwas_catalog,
        gwas_catalog_ids,
        my_snp_data,
        my_snp_run_date,
        chrom,
        start,
        end,
        snp_finemap_signals,
        str_finemap_signals,
        finemap_regions,
        conditioned_strs,
        conditioned_isnps,
        *,
        return_figure = False):
    ext = outfname.split('.')[-1]

    if not binary:
        stat = 'mean'
    else:
        stat = 'fraction'

    conditioning = conditioned_strs or conditioned_isnps
    if conditioned_isnps:
        print("Can't handle conditioned_isnps at the moment")
        exit()

    if return_figure:
        assert outfname[0] == '.'
        assert ext != 'html'

    plot_my_str_data = my_str_data is not None
    plot_my_snp_data = my_snp_data is not None
    plot_plink_snp_data = plink_snp_data is not None
    plot_gwas_catalog = gwas_catalog is not None
    if plot_my_str_data:
        plot_peaks = 'is_peak' in my_str_data.dtype.names
    else:
        plot_peaks = False

    assert bool(snp_finemap_signals) + bool(str_finemap_signals) + bool(finemap_regions) in {0, 3}

    print(f"Plotting phenotype {phenotype} ... ", end='', flush=True)
    start_time = time.time()
    # reformat data for plotting

    cols_to_skip = {
        'coeff_phenotype',
        'coeff_intercept',
        'R^2',
        'total_hardcall_genotypes',
        'subset_total_hardcall_genotypes'
    }

    sources = []
    source_dicts = []

    if not conditioning and not start:
        start_p_val_cap = 30
    else:
        start_p_val_cap = None

    if plot_my_str_data:
        my_str_source, my_str_sources = create_source_dict(
            my_str_data,
            chrom,
            cols_to_skip = cols_to_skip,
            chrs_to_var_signals = str_finemap_signals,
            start_p_val_cap = start_p_val_cap
        )
        sources.append(my_str_source)
        source_dicts.append(my_str_sources)
        for source in (my_str_source, *my_str_sources.values()):
            source.data['size'] = np.full(source.data['pos'].shape, 6)
            if plot_peaks:
                source.data['size'][source.data['is_peak']] = 15

    if plot_my_snp_data:
        my_snp_source, my_snp_sources = create_source_dict(
            my_snp_data,
            chrom,
            cols_to_skip = cols_to_skip,
            start_p_val_cap = start_p_val_cap
        )
        sources.append(my_snp_source)
        source_dicts.append(my_snp_sources)

    if plot_plink_snp_data:
        plink_snp_data = numpy.lib.recfunctions.merge_arrays(
            (
                plink_snp_data,
                np.rec.fromarrays((np.char.add(np.char.add(
                    plink_snp_data['ref'], ','), plink_snp_data['alt']
                ),), names=['alleles'])
            ), flatten=True
        )
        plink_snp_source, plink_snp_sources = create_source_dict(
            plink_snp_data,
            chrom,
            chrs_to_var_signals = snp_finemap_signals
        )
        sources.append(plink_snp_source)
        source_dicts.append(plink_snp_sources)
        for source in (plink_snp_source, *plink_snp_sources.values()):
            source.data['size'] = np.full(source.data['pos'].shape, 5)
            if plot_peaks:
                source.data['size'][source.data['is_peak']] = 15

    if plot_gwas_catalog:
        catalog_source, catalog_sources = create_source_dict(np.rec.fromarrays((
            gwas_catalog[:, 0],
            gwas_catalog[:, 1],
            gwas_catalog[:, 2],
            gwas_catalog_ids
        ), names=['chr', 'pos', 'p_val', 'rsids']),
        chrom)
        sources.append(catalog_source)
        source_dicts.append(catalog_sources)

    if plot_my_str_data and ext == 'html':
        locus_plot = bokeh.plotting.figure(
            width=400,
            height=400,
            title='Click on an STR in the Manhattan plot!                                      ',
            x_axis_label='Sum of allele lengths (repeat copies)',
            y_axis_label=f'Mean {phenotype} ({unit})',
            tools='save',
            y_range=(-2, 5)
        )
        subtext = "Phenotype values are unadjusted for covariates"
        if conditioning:
            subtext += " or genotypes that were conditioned on"
        locus_plot.add_layout(bokeh.models.Title(
            text=subtext,
            align='right'
        ), 'below')

        locus_source = bokeh.models.ColumnDataSource({
            'x': [1, 2],
            'y': [1, 2],
            'CI5e_2_lower': ['0', '1'],
            'CI5e_2_upper': ['3', '4'],
            'CI5e_8_lower': ['-1', '-2'],
            'CI5e_8_upper': ['4', '4.5']
        })

        locus_plot.varea('x', 'CI5e_2_upper', 'CI5e_8_upper', source=locus_source, color="red", alpha=0.2)
        locus_plot.varea('x', 'CI5e_2_upper', 'CI5e_2_lower', source=locus_source, color="red", alpha=0.4)
        locus_plot.varea('x', 'CI5e_2_lower', 'CI5e_8_lower', source=locus_source, color="red", alpha=0.2)
        locus_plot.line('x', 'y', source=locus_source, line_width=2, color="black")
        locus_plot.circle('x', 'y', source=locus_source, color="black", size=6)
        
        locus_wheel_zoom = bokeh.models.tools.WheelZoomTool(dimensions="height")
        locus_plot.add_tools(locus_wheel_zoom)
        locus_plot.toolbar.active_scroll = locus_wheel_zoom

        locus_pan = bokeh.models.tools.PanTool(dimensions="height")
        locus_plot.add_tools(locus_pan)
        locus_plot.toolbar.active_drag = locus_pan

    # set up drawing canvas
    if chrom is None:
        x_axis_label = 'Position (bp)'
    else:
        x_axis_label = f'Position (bp chr{chrom})'

    if ext == 'html':
        output_backend = 'webgl'
    elif ext == 'svg':
        output_backend = 'svg'
    else:
        assert ext == 'png'
        output_backend = 'canvas'

    if start_p_val_cap is not None:
        start_height_cap = start_p_val_cap
    else:
        if plot_plink_snp_data:
            plink_data = plink_snp_source.data
            start_height_cap = np.max(
                plink_data['p_val'][(start <= plink_data['pos']) & (plink_data['pos'] <= end)]
            )

        if plot_my_str_data:
            str_data = my_str_source.data
            included_strs = str_data['p_val'][(start <= str_data['pos']) & (str_data['pos'] <= end)]
            if included_strs.shape[0] > 0:
                start_height_cap = max(start_height_cap, np.max(included_strs))

    manhattan_plot = bokeh.plotting.figure(
        width=1200,
        height=900,
        title=(phenotype.capitalize() + ' Manhattan Plot'),
        x_axis_label=x_axis_label,
        y_axis_label='-log10(p-value)',
        tools='xzoom_in,xzoom_out,save',
        y_range=(-0.025*start_height_cap, start_height_cap*1.025),
        output_backend=output_backend
    )
    if start:
        x_width = end - start
        manhattan_plot.x_range = bokeh.models.Range1d(start-0.025*x_width, end+0.025*x_width)

    # add custom tools
    wheel_zoom = bokeh.models.tools.WheelZoomTool(dimensions="width")
    manhattan_plot.add_tools(wheel_zoom)
    manhattan_plot.toolbar.active_scroll = wheel_zoom

    pan = bokeh.models.tools.PanTool(dimensions="width")
    manhattan_plot.add_tools(pan)
    manhattan_plot.toolbar.active_drag = pan

    # make sure only one tooltip displays, not all moused over tooltips
    # from https://stackoverflow.com/a/64544102/2966505
    # and https://stackoverflow.com/a/707580/2966505
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
    manhattan_plot.js_on_event(bokeh.events.MouseEnter, one_tooltip_callback)

    if not start:
        start = 0
        if not chrom:
            end = chr_lens[1]
        else:
            end = chr_lens[chrom]
    line_source = bokeh.models.ColumnDataSource(dict(
        x=[start, end],
        y=[-np.log10(5e-8)]*2,
    ))
    manhattan_plot.line(
        x='x',
        y='y',
        source=line_source,
        line_width=3,
        color='red',
        legend_label='GWAS p-value threshold'
    )

    # colors from here: https://bconnelly.net/posts/creating_colorblind-friendly_figures/
    # RGB255
    # 0,0,0
    # 230, 159, 0
    # 86, 180, 233
    # 0, 158, 115
    # 240, 228, 66 - yellow, draws attention and has some problems
    # 0, 114, 178
    # 213, 94, 0
    # 204, 121, 167
    if finemap_regions:
        cmap_field_name = 'FINEMAP_pcausal'
    else:
        cmap_field_name = 'pos' # arbitrary existant field

    if plot_plink_snp_data:
        plink_snp_color = bokeh.colors.RGB(0, 114, 178)
        plink_snp_cmap = bokeh.transform.linear_cmap(
            field_name = cmap_field_name,
            low = 0,
            high = 0.2,
            palette = colorcet.kr,
            nan_color = 'purple'
        )
        plink_snp_manhattan = manhattan_plot.circle(
            'pos',
            'display_p_val',
            source=plink_snp_source,
            legend_label='SNPs Plink',
            color=plink_snp_color,
            size='size',
            muted_alpha=0.1
        )
    if plot_my_str_data:
        my_str_color = bokeh.colors.RGB(204, 121, 167)
        my_str_cmap = bokeh.transform.linear_cmap(
            field_name = cmap_field_name,
            low = 0,
            high = 0.2,
            palette = colorcet.kr,
            nan_color = 'purple'
        )
        my_str_manhattan = manhattan_plot.square_pin(
            'pos',
            'display_p_val',
            source=my_str_source,
            legend_label='STRs my code',
            color=my_str_color,
            size='size',
            muted_alpha=0.1
        )
    if plot_my_snp_data:
        my_snp_manhattan = manhattan_plot.circle(
            'pos',
            'display_p_val',
            source=my_snp_source,
            legend_label='SNPs my code',
            color=(204, 121, 167),
            size=5,
            muted_alpha=0.1
        )
    if plot_gwas_catalog:
        catalog_manhattan = manhattan_plot.square(
            'pos',
            'display_p_val',
            source=catalog_source,
            legend_label='GWAS Catalog Hits',
            color=(213, 94, 0),
            size=7,
            muted_alpha=0.1
        )
    if conditioned_isnps:
        pass #TODO
    if conditioned_strs:
        manhattan_plot.square_pin(
            [int(STR) for STR in conditioned_strs],
            [0]*len(conditioned_strs),
            legend_label='Conditioned-on STRs',
            color=(0, 158, 115),
            size=20
        )

    # add hover tooltips for each manhattan plot
    # see https://stackoverflow.com/questions/49282078/multiple-hovertools-for-different-lines-bokeh
    hover_tools = []
    if plot_my_str_data and ext == 'html':
        my_str_hover = bokeh.models.tools.HoverTool(renderers=[my_str_manhattan])
        my_str_hover.tooltips = [
            ('var type', 'STR'),
            ('alleles:', '@alleles'),
            ('pos', '@pos'),
            ('-log10(p_val) my code', '@p_val')
        ]
        if plot_peaks:
            my_str_hover.tooltips.append(
                ('is there a nearby SNP in this peak?', '@tagged_by_other_variant_type')
            )
        if binary:
            my_str_hover.tooltips.append(
                ('Using Firth penalization during regression?', '@firth')
            )
        if conditioning:
            my_str_hover.tooltips.append(
                ('-log10(unconditioned_p_val) my code', '@unconditioned_p')
            )
        if str_finemap_signals:
            my_str_hover.tooltips.append(('FINEMAP_pcausal', '@FINEMAP_pcausal{safe}'))
        str_stats_start_idx = list(my_str_data.dtype.names).index(
            f'{stat}_{phenotype}_per_single_dosage'
        )
        for detail_name in my_str_data.dtype.names[7:str_stats_start_idx]:
            if detail_name in cols_to_skip:
                continue
            my_str_hover.tooltips.append((detail_name, f'@{detail_name}' '{safe}'))
        manhattan_plot.add_tools(my_str_hover)
        hover_tools.append(my_str_hover)

    if plot_my_snp_data and ext == 'html':
        my_snp_hover = bokeh.models.tools.HoverTool(renderers=[my_snp_manhattan])
        my_snp_hover.tooltips = [
            ('var type', 'STR'),
            ('alleles:', '@alleles'),
            ('pos', '@pos'),
            ('-log10(p_val) my code', '@p_val')
        ]
        snp_stats_start_idx = list(my_snp_data.dtype.names).index(
            f'{stat}_{phenotype}_per_single_dosage'
        )
        for detail_name in my_snp_data.dtype.names[7:snp_stats_start_idx]:
            if detail_name in cols_to_skip:
                continue
            my_snp_hover.tooltips.append((detail_name, f'@{detail_name}' '{safe}'))
        manhattan_plot.add_tools(my_snp_hover)
        hover_tools.append(my_snp_hover)

    if plot_plink_snp_data and ext == 'html':
        plink_snp_hover = bokeh.models.tools.HoverTool(renderers=[plink_snp_manhattan])
        plink_snp_hover.tooltips = [
            ('var type', 'SNP'),
            ('alleles:', '@alleles'),
            ('pos', '@pos'),
            ('ID', '@id'),
            ('-log10(p_val) Plink', '@p_val')
        ]
        if binary == 'logistic':
            plink_snp_hover.tooltips.append(
                ('alt_case_count', '@alt_case_count')
            )
            plink_snp_hover.tooltips.append(
                ('alt_control_count', '@alt_control_count')
            )
            plink_snp_hover.tooltips.append(
                ('Using Firth penalization during regression?', '@firth')
            )
        if plot_peaks:
            plink_snp_hover.tooltips.append(
                ('is there a nearby STR in this peak?', '@tagged_by_other_variant_type')
            )
        if conditioning:
            plink_snp_hover.tooltips.append(
                ('-log10(unconditioned_p_val) Plink', '@unconditioned_p')
            )
            plink_snp_hover.tooltips.append(
                ('plink error code', '@error')
            )
        plink_snp_hover.tooltips.extend([
            ('Minor allele frequency', '@maf'),
            ('Imputation INFO', '@info')
        ])
        if snp_finemap_signals:
            plink_snp_hover.tooltips.append(('FINEMAP_pcausal', '@FINEMAP_pcausal' '{safe}'))
        manhattan_plot.add_tools(plink_snp_hover)
        hover_tools.append(plink_snp_hover)

    if plot_gwas_catalog:
        catalog_hover = bokeh.models.tools.HoverTool(renderers=[catalog_manhattan])
        catalog_hover.tooltips = [
            ('var type', 'GWAS Catalog SNP'),
            ('pos', '@pos'),
            ('rsid:', '@rsids'),
            ('-log10(p_val)', '@p_val')
        ]
        manhattan_plot.add_tools(catalog_hover)
        hover_tools.append(catalog_hover)

    manhattan_plot.toolbar.active_inspect = hover_tools

    if plot_my_str_data and ext == 'html':
        tap = bokeh.models.tools.TapTool()
        manhattan_plot.add_tools(tap)
        manhattan_plot.toolbar.active_tap = tap

        plot_tap_callback = bokeh.models.CustomJS(
            args=dict(
                locus_plot=locus_plot,
                y_range=locus_plot.y_range,
                x_range=locus_plot.x_range,
                my_str_source=my_str_source,
                locus_source=locus_source
            ),
            code = f"""
                if (my_str_source.selected.indices.length == 0) {{
                    return;
                }}
                var idx = my_str_source.selected.indices[0];
                var data = my_str_source.data;
                var locus_dict = JSON.parse(data['{stat}_{phenotype}_per_single_dosage'][idx]);
                var CI5e_2_str = data['CI5e_2SingleDosagePhenotype'][idx].replaceAll('(', '[').replaceAll(')', ']').replaceAll('nan', '"NaN"');
                var CI5e_2_dict = JSON.parse(CI5e_2_str);
                var CI5e_8_str = data['CI5e_8SingleDosagePhenotype'][idx].replaceAll('(', '[').replaceAll(')', ']').replaceAll('nan', '"NaN"');
                var CI5e_8_dict = JSON.parse(CI5e_8_str);
                var gts = [];
                var str_gts = [];
                var phenotypes = [];
                var CI5e_2_lower = [];
                var CI5e_2_upper = [];
                var CI5e_8_lower = []; 
                var CI5e_8_upper = []; 
                Object.entries(locus_dict).forEach(([key, value]) => {{
                    // if (value != 'NaN') {{
                        str_gts.push(key);
                        gts.push(parseFloat(key));
                        phenotypes.push(parseFloat(value));
                    // }}
                }});
                var y_max = Math.max.apply(Math, phenotypes);
                var y_min = Math.min.apply(Math, phenotypes);
                y_range.start = y_min - 0.05*(y_max - y_min);
                y_range.end = y_max + 0.05*(y_max - y_min);
                y_range.change.emit();
                var x_max = Math.max.apply(Math, gts);
                var x_min = Math.min.apply(Math, gts);
                x_range.start = x_min - 0.05*(x_max - x_min);
                x_range.end = x_max + 0.05*(x_max - x_min);
                x_range.change.emit();
                for (var i=0; i<gts.length; i++) {{
                    CI5e_2_lower.push(CI5e_2_dict[str_gts[i]][0]); 
                    CI5e_2_upper.push(CI5e_2_dict[str_gts[i]][1]); 
                    CI5e_8_lower.push(CI5e_8_dict[str_gts[i]][0]); 
                    CI5e_8_upper.push(CI5e_8_dict[str_gts[i]][1]); 
                }}
                
                console.log(
                    gts,
                    phenotypes,
                    CI5e_2_lower, 
                    CI5e_2_upper, 
                    CI5e_8_lower, 
                    CI5e_8_upper
                )
                var new_dict = {{
                    'x': gts,
                    'y': phenotypes,
                    'CI5e_2_lower': CI5e_2_lower,
                    'CI5e_2_upper': CI5e_2_upper,
                    'CI5e_8_lower': CI5e_8_lower,
                    'CI5e_8_upper': CI5e_8_upper
                }};
                locus_source.data = new_dict; 
                locus_source.change.emit();

                locus_plot.title.text = 'STR ' + data['chr'][idx].toString() + ':' + data['pos'][idx].toString() + ' Repeat Unit: ' + data['motif'][idx] + ' -log10(Association p-val): ' + data['p_val'][idx].toFixed(2);
            """
        )
        my_str_source.selected.js_on_change('indices', plot_tap_callback)

    height_slider = bokeh.models.Slider(
        start = 8,
        end=max_p_val,
        value=start_height_cap,
        step=1,
        title="p-value cap",
        sizing_mode = 'stretch_width'
    )
    for source in sources:
        height_callback = bokeh.models.CustomJS(
            args=dict(source=source, height_slider=height_slider),
            code="""
                const data = source.data;
                const display_p_val = data['display_p_val'];
                const p_val= data['p_val'];
                const new_max = height_slider.value;
                for (var i = 0; i < display_p_val.length; i++) {
                    display_p_val[i] = Math.min(p_val[i], new_max);
                }
                source.change.emit();
            """
        )
        height_slider.js_on_change('value', height_callback)
    height_callback = bokeh.models.CustomJS(
        args=dict(height_slider=height_slider, y_range=manhattan_plot.y_range),
        code="""
            const new_max = height_slider.value;
            y_range.start = -0.025*new_max
            y_range.end = new_max*1.025;
            y_range.change.emit();
        """
    )
    height_slider.js_on_change('value', height_callback)

    if chrom is None:
        chrom_select = bokeh.models.Select(
            title="Chromosome:",
            options=[str(num) for num in range(1, 23)]
        )
        for source, source_dict in zip(sources, source_dicts):
            chrom_callback = bokeh.models.CustomJS(
                args=dict(
                    source=source,
                    source_dict=source_dict,
                    height_slider = height_slider,
                    line_source = line_source,
                    chr_lens = chr_lens,
                    range=manhattan_plot.x_range
                ),
                code = """
                    source.data = source_dict[this.value].data;
                    const data = source.data;
                    const display_p_val = data['display_p_val'];
                    const p_val= data['p_val'];
                    const new_max = height_slider.value;
                    for (var i = 0; i < display_p_val.length; i++) {
                        display_p_val[i] = Math.min(p_val[i], new_max);
                    }
                    source.change.emit();
                    line_source.data['x'] = [0, chr_lens[this.value-1]];
                    line_source.change.emit();
                    range.start = 0;
                    range.end = chr_lens[this.value-1];
                """
            )
            chrom_select.js_on_change('value', chrom_callback)

    if snp_finemap_signals:
        finemap_toggle = bokeh.models.Toggle(label="Color by FINEMAP causal posterior probability")
        finemap_toggle.js_on_click(bokeh.models.CustomJS(
            args=dict(
                sources = [my_str_source, plink_snp_source],
                manhattans = [my_str_manhattan, plink_snp_manhattan],
                colors = [my_str_color, plink_snp_color],
                cmaps = [my_str_cmap, plink_snp_cmap]
            ),
            code = """
                if (this.active) { for (var i=0; i<cmaps.length; i++) {
                    manhattans[i].glyph.fill_color = cmaps[i];
                    manhattans[i].glyph.line_color = cmaps[i];
                    sources[i].change.emit();
                }} else { for (var i=0; i<cmaps.length; i++) {
                    manhattans[i].glyph.fill_color = colors[i];
                    manhattans[i].glyph.line_color = colors[i];
                    sources[i].change.emit();
                }}
            """
        ))

        use_chrom = chrom
        if use_chrom is None:
            use_chrom = 1
        finemap_region_select = bokeh.models.Select(
            title="Select finemap region:",
            options=finemap_regions[use_chrom]
        )
        if chrom is None:
            finemap_region_select.js_on_change('value', bokeh.models.CustomJS(
                args = dict(
                    x_range=manhattan_plot.x_range,
                    chr_lens=chr_lens,
                    chrom_select=chrom_select
                ),
                code = """
                    var chrom_max = chr_lens[parseInt(chrom_select.value) - 1];
                    var pieces = this.value.split('_');
                    x_range.start = Math.max(parseInt(pieces[0]) - 50000, 0);
                    x_range.end = Math.min(parseInt(pieces[1]) + 50000, chrom_max);
                    x_range.change.emit();
                """
            ))
            chrom_select.js_on_change('value', bokeh.models.CustomJS(
                args=dict(
                    region_select=finemap_region_select,
                    regions=finemap_regions
                ),
                code = """
                    region_select.value = '';
                    region_select.options = regions[parseInt(this.value)];
                    region_select.change.emit();
                """
            ))
        else:
            finemap_region_select.js_on_change('value', bokeh.models.CustomJS(
                args = dict(
                    x_range=manhattan_plot.x_range,
                    chrom_max=chr_lens[chrom - 1]
                ),
                code = """
                    var pieces = this.value.split('_');
                    x_range.start = Math.max(parseInt(pieces[0]) - 50000, 0);
                    x_range.end = Math.min(parseInt(pieces[1]) + 50000, chrom_max);
                    x_range.change.emit();
                """
            ))


    if chrom is None:
        chrom_layout_row = [height_slider, chrom_select]
    else:
        chrom_layout_row = [height_slider]

    if snp_finemap_signals:
        chrom_layout_row.append(finemap_region_select)
        chrom_layout_row.append(finemap_toggle)

    if ext == 'html':
        if plot_my_str_data:
            layout = bokeh.layouts.grid([
                locus_plot,
                chrom_layout_row,
                manhattan_plot
            ])
        else:
            layout = bokeh.layouts.grid([
                chrom_layout_row,
                manhattan_plot
            ])

    if my_str_run_date:
        manhattan_plot.add_layout(bokeh.models.Title(
            text=f"My STR code run date: {my_str_run_date}",
            align='right'
        ), 'below')
    if my_snp_run_date:
        manhattan_plot.add_layout(bokeh.models.Title(
            text=f"My SNP code run date: {my_snp_run_date}",
            align='right'
        ), 'below')
    if plink_snp_run_date:
        manhattan_plot.add_layout(bokeh.models.Title(
            text=f"Plink SNP run date: {plink_snp_run_date}",
            align='right'
        ), 'below')

    manhattan_plot.legend.click_policy="mute"

    if return_figure:
        return manhattan_plot

    if ext == 'html':
        html = bokeh.embed.file_html(layout, bokeh.resources.CDN, f'Manhattan plot {phenotype}')
        with open(outfname, 'w') as outfile:
            outfile.write(html)
    elif ext == 'svg':
        bokeh.io.export_svg(manhattan_plot, filename=outfname)
    else:
        assert ext == 'png'
        manhattan_plot.background_fill_color = None
        manhattan_plot.border_fill_color = None
        bokeh.io.export_png(manhattan_plot, filename=outfname)
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

my_results_rename = {
    0: 'chr',
    3: 'filtered',
    4: 'p_val',
    5: 'coeff_phenotype'
}
my_str_results_rename = {
    -5: 'CI5e_2SingleDosagePhenotype',
    -4: 'CI5e_8SingleDosagePhenotype',
    -2: 'CI5e_2PairedDosagePhenotype',
    -1: 'CI5e_8PairedDosagePhenotype'
}

def load_data(phenotype, binary, condition):
    if not condition:
        return (
            load_my_str_results(phenotype, binary, None),
            load_plink_results(phenotype, binary, None),
            *load_gwas_catalog(phenotype)
        )
    else:
        return (
            load_my_str_results(phenotype, binary, condition),
            load_plink_results(phenotype, binary, condition),
            None,
            None
        )

def load_my_str_results(phenotype, binary, condition):
    print(f"Loading my STR results for {phenotype} ... ", end='', flush=True)
    if binary:
        runtype_suffix = '_' + binary
    else:
        runtype_suffix = ''
    start_time = time.time()
    unconditioned_results_fname = f'{ukb}/association/plots/input/{phenotype}/my_str{runtype_suffix}_results.tab'
    unconditioned_results = pd.read_csv(
        unconditioned_results_fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        dtype=utils.get_dtypes(unconditioned_results_fname, {'locus_filtered': str})
    )

    if not condition:
        results = unconditioned_results
    else:
        conditional_results_fname = \
            f'{ukb}/association/results/{phenotype}/my_str{runtype_suffix}_conditional/{condition}.tab'

        results = pd.read_csv(
            conditional_results_fname,
            header=0,
            delimiter='\t',
            encoding='UTF-8',
            dtype=utils.get_dtypes(conditional_results_fname, {'locus_filtered': str})
        )

        unconditioned_results[f'p_{phenotype}'] = np.maximum(unconditioned_results[f'p_{phenotype}'], 1 / 10**max_p_val)
        unconditioned_results[f'p_{phenotype}'] = -np.log10(unconditioned_results[f'p_{phenotype}'])
        unconditioned_results.rename(
            columns = {f'p_{phenotype}': 'unconditioned_p'}, inplace=True
        )
        unconditioned_results = unconditioned_results[['chrom', 'pos', 'unconditioned_p']]

        results = results.merge(
            unconditioned_results,
            on=['chrom', 'pos'],
            how= 'inner'
        ) # subsets to only those which passed the p-val threshold in the unconditioned run

    if binary == 'logistic':
        results.rename(columns={'firth?': 'firth'}, inplace=True)

    rename_dict = {}
    for idx, name in my_results_rename.items():
        rename_dict[results.columns[idx]] = name
    for idx, name in my_str_results_rename.items():
        if condition:
            idx -=1
        rename_dict[results.columns[idx]] = name
    for colname in ('total_per_allele_dosages', 'total_hardcall_alleles',
                'subset_total_per_allele_dosages', 'subset_total_hardcall_alleles',
                'subset_allele_dosage_r2'):
        # convert allele lens from strings to floats, in addition round allele lens and values, but not NaN values
        new_col = np.array(list(map(
            lambda dict_str: {round(float(allele_len), 2): (round(val, 2) if val != 'NaN' else val) for allele_len, val in ast.literal_eval(dict_str).items()},
            results[colname]
        )))
        # convert allele_lens to ints if they are close enough
        new_col = np.array(list(map(
            lambda d: str({(int(key) if key == int(key) else key) : val for key, val in d.items()}),
            new_col
        )))
        results[colname] = new_col
    results.rename(columns=rename_dict, inplace=True)
    results = utils.df_to_recarray(results)
    results['p_val'] = np.maximum(results['p_val'], 1 / 10**max_p_val)
    results['p_val'] = -np.log10(results['p_val'])
    if condition:
        for STR in get_conditioned_strs(condition):
            results = results[results['pos'] != STR]
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return results

def load_my_snp_results(phenotype, binary, chrom, start, end):
    print(f"Loading my {phenotype} snp results chrom={chrom} start={start} end={end} ... ", end='', flush=True)
    start_time = time.time()

    snp_results_fname = f'{ukb}/association/plots/input/{phenotype}/my_imputed_snp'
    if binary:
        snp_results_fname += '_' + binary
    snp_results_fname += f'_chr{chrom}'
    if start:
        snp_results_fname += f'_{start}_{end}'
    snp_results_fname += '_results.tab'

    type_overrides = {
        'locus_filtered': str,
        f'p_{phenotype}': float,
        '0.05_significance_CI': str,
        '5e-8_significance_CI': str
    }
    if not binary:
        type_overrides.update({f'mean_{phenotype}_per_single_dosage': str})
    else:
        type_overrides.update({f'fraction_{phenotype}_per_single_dosage': str})

    if binary == 'logistic':
        type_overrides.update({'unused_col': str, 'firth?': str})

    my_snp_results = utils.df_to_recarray(pd.read_csv(
        snp_results_fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        dtype=utils.get_dtypes(snp_results_fname, type_overrides)
    ))

    names = list(my_snp_results.dtype.names)
    for idx, name in my_results_rename.items():
        names[idx] = name
    my_snp_results.dtype.names = names

    my_snp_results['p_val'] = -np.log10(my_snp_results['p_val'])

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return my_snp_results

def load_plink_results(phenotype, binary, condition, fname=None):
    # TODO remove conditioned snps
    # Load plink SNP results
    print(f"Loading plink SNP results for {phenotype} ... ", end='', flush=True)
    start_time = time.time()

    assert bool(condition) + bool(fname) <= 1

    if fname is None:
        if binary:
            runtype_suffix = '_' + binary
        else:
            runtype_suffix = ''
        fname = f'{ukb}/association/plots/input/{phenotype}/plink_snp{runtype_suffix}_results_with_mfi.npy'

    unconditioned_results = pd.DataFrame.from_records(np.load(fname))
    if not condition:
        results = unconditioned_results
    else:
        results = pd.DataFrame.from_records(np.load(
            f'{ukb}/association/plots/input/{phenotype}/plink_snp{runtype_suffix}_conditional_{condition}_results_with_mfi.npy'
        ))

        unconditioned_results['p_val'] = np.maximum(unconditioned_results['p_val'], 1 / 10**max_p_val)
        unconditioned_results['p_val'] = -np.log10(unconditioned_results['p_val'])
        unconditioned_results.rename(
            columns={'p_val': 'unconditioned_p'},
            inplace=True
        )
        unconditioned_results = unconditioned_results[['chr', 'pos', 'unconditioned_p']]

        results = results.merge(
            unconditioned_results,
            on=['chr', 'pos'],
            how = 'inner'
        ) # subsets to only those which passed the p-val threshold in the unconditioned run

    if binary == 'logistic':
        results.rename(columns={'firth?': 'firth'}, inplace=True)

    results = utils.df_to_recarray(results)

    results = results[results['error'] != 'CONST_OMITTED_ALLELE']
    if binary == 'logistic':
        # in theory could keep unfinished error codes and just note them,
        # but easier to ignore
        results = results[
            (results['error'] != 'FIRTH_CONVERGE_FAIL') &
            (results['error'] != 'UNFINISHED')
        ]
    results['p_val'] = np.maximum(results['p_val'], 1 / 10**max_p_val)
    results['p_val'] = -np.log10(results['p_val'])

    # we've already filtered all the spots that had errors in the unconditional run
    # having a VIF_TOO_HIGH or CORR_TOO_HIGH only in the conditional run just means that
    # SNP is extremely correlated with the conditioning variants, which means
    # its p-value should be very small, so this isn't an issue.
    if not condition:
        if not np.all(results['error'] == '.'):
            print(np.unique(results['error']))
            assert False
    else:
        assert np.all(
            (results['error'] == '.') |
            (results['error'] == 'VIF_TOO_HIGH') |
            (results['error'] == 'CORR_TOO_HIGH')
        )
    # rename for readability
    results['error'][results['error'] == '.'] = 'none'
    results['p_val'][results['error'] == 'VIF_TOO_HIGH'] = 0

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return results

def load_gwas_catalog(phenotype):
    if phenotype not in {'height', 'total_bilirubin'}:
        return None, None
    # Load the NHGRI-EBI GWAS catalog
    print("Loading GWAS catalog results ... ", end='', flush=True)
    start_time = time.time()
    catalog = np.loadtxt(
        f'{ukb}/misc_data/snp_summary_stats/catalog/catalog_hg19.tsv',
        usecols=[11, 12, 27],
        delimiter='\t',
        skiprows=1,
        dtype=object
    ) # chrom, pos, pvalue
    # omit results that are not mapped to a recognizable chromosome
    catalog_names = np.loadtxt(
        f'{ukb}/misc_data/snp_summary_stats/catalog/catalog_hg19.tsv',
        usecols=[7, 21, 34],
        delimiter='\t',
        skiprows=1,
        dtype=object
    ).astype('U') # cols 7, 34 describe the 'height', 21 is the snp rsid
    catalog_names = np.char.lower(catalog_names)
    filter_weird_chroms = np.isin(catalog[:, 0],
                                  list(str(chrom) for chrom in range(1, 23)))
    catalog_names = catalog_names[filter_weird_chroms, :]
    catalog = catalog[filter_weird_chroms, :]
    catalog = catalog.astype(float)

    known_assocs = {}
    known_assoc_ids = {}
    height_rows = np.logical_and(
        catalog_names[:, 2] == 'body height',
        catalog_names[:, 0] != "pericardial adipose tissue adjusted for height and weight"
    )
    known_assocs['height'] = catalog[height_rows, :]
    known_assoc_ids['height'] = catalog_names[height_rows, 1]
    known_assocs['total_bilirubin'] = \
            catalog[catalog_names[:, 2] == 'bilirubin measurement', :]
    known_assoc_ids['total_bilirubin'] = \
            catalog_names[catalog_names[:, 2] == 'bilirubin measurement', 1]

    for assoc_phen in known_assocs:
        zero_idx = known_assocs[assoc_phen][:, 2] == 0
        known_assocs[assoc_phen][~zero_idx, 2] = -np.log10(known_assocs[assoc_phen][~zero_idx, 2])
        known_assocs[assoc_phen][zero_idx, 2] = 50.12345 #choose an arbitrary number

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

    return known_assocs[phenotype], known_assoc_ids[phenotype]

def load_finemap_signals(finemap_signals):
    # finemap signals must be sorted
    print("Loading FINEMAP causality results ... ", end='', flush=True)
    start_time = time.time()
    chrs_to_snp_signals = {}
    chrs_to_str_signals = {}
    chrs_to_regions = {}
    for i in range(1, 23):
        chrs_to_snp_signals[i] = []
        chrs_to_str_signals[i] = []
        chrs_to_regions[i] = []
    for signal in finemap_signals:
        region = signal.split('/')[-2]
        chrom, start, end = (int(val) for val in region.split('_'))
        chrs_to_regions[chrom].append(f'{start}_{end}')
        with open(signal) as per_var_output:
            next(per_var_output)
            for line in per_var_output:
                _id, pos, pcausal = np.array(line.split())[[1, 3, 10]]
                pos = int(pos)
                split = _id.split('_')
                assert int(split[1]) == pos
                if pcausal == 'NA':
                    pcausal = np.nan
                else:
                    pcausal = float(pcausal)
                if _id[:4] == 'STR_':
                    chrs_to_str_signals[chrom].append((pos, pcausal))
                elif _id[:4] == 'SNP_':
                    chrs_to_snp_signals[chrom].append(
                        (pos, split[2], split[3], pcausal)
                    )
                else:
                    raise ValueError(f'Found uninterpretable id {_id}')

    chrs_to_snp_signals = {
        key: (
            pd.DataFrame(np.stack(val), columns=('pos', 'ref', 'alt', 'pcausal'))
            if len(val) > 0
            else pd.DataFrame(      [], columns=('pos', 'ref', 'alt', 'pcausal'))
        )
        for key, val in chrs_to_snp_signals.items()
    }
    chrs_to_str_signals = {
        key: (
            pd.DataFrame(np.stack(val), columns=('pos', 'pcausal'))
            if len(val) > 0
            else pd.DataFrame(      [], columns=('pos', 'pcausal'))
        )
        for key, val in chrs_to_str_signals.items()
    }

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return (chrs_to_snp_signals, chrs_to_str_signals, chrs_to_regions)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('ext', help='file extension', choices=['png', 'svg', 'html'])
    parser.add_argument("--my-plink-comparison", action='store_true', default=False)
    parser.add_argument("--only-plink", action='store_true', default=False)
    parser.add_argument("--only-mine", action='store_true', default=False)
    parser.add_argument("--chrom", type=int)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--condition")
    parser.add_argument('--finemap-signals', action='store_true', default=False)
    parser.add_argument('--binary', default=False, choices={'linear', 'logistic'})
    parser.add_argument('--plink-snp-fname')
    parser.add_argument('--peaks-spacing', type=int)
    parser.add_argument('--peaks-thresh', type=str)
    parser.add_argument('--overview', default=False, action='store_true')
    args = parser.parse_args()

    phenotype = args.phenotype
    ext = args.ext
    binary = args.binary
    condition = args.condition

    if not args.binary:
        run_suffix = ''
    else:
        run_suffix = '.' + binary

    assert args.finemap_signals + args.my_plink_comparison + bool(condition) + bool(args.peaks_spacing) <= 1

    assert bool(args.peaks_spacing) == bool(args.peaks_thresh)
    if args.overview:
        assert bool(args.peaks_spacing)
        assert not bool(args.chrom)
        assert ext == 'png'

    assert bool(args.start) == bool(args.end)
    if bool(args.start):
        assert bool(args.chrom)

    if bool(args.condition):
        # will parse chrom, start, end from the condition string itself
        assert not bool(args.chrom)
        assert ext == 'png'
    else:
        conditioned_isnps = conditioned_strs = None

    if args.my_plink_comparison:
        assert bool(args.chrom)
    if bool(args.plink_snp_fname) or args.only_mine or args.only_plink:
        assert args.my_plink_comparison
    assert args.only_mine + args.only_plink <= 1

    if not args.finemap_signals:
        snp_finemap_signals = None
        str_finemap_signals = None
        finemap_regions = None

    if not args.my_plink_comparison and not args.condition:
        with open(f'{ukb}/association/results/{phenotype}/my_str{run_suffix.replace(".", "_")}/README.txt') as README:
            date_line = next(README)
            my_str_run_date = date_line.split(' ')[2]
    else:
        my_str_run_date = None

    if not bool(args.my_plink_comparison):
        my_str_results, plink_snp_results, gwas_catalog, gwas_catalog_ids = load_data(
            phenotype, binary, condition
        )
        my_snp_results = None
        my_snp_run_date = None
        if not condition and not args.chrom:
            chrom = None
            start = None
            end = None
            outfname = f'{ukb}/association/plots/{phenotype}{run_suffix}.manhattan.{ext}'
        elif args.chrom:
            chrom = args.chrom
            if args.start:
                start = args.start
                end = args.end
                outfname = f'{ukb}/association/plots/{phenotype}{run_suffix}.manhattan.chr{chrom}_{start}_{end}.{ext}'
            else:
                start = None
                end = None
                outfname = f'{ukb}/association/plots/{phenotype}{run_suffix}.manhattan.chr{chrom}.{ext}'
        else:
            chrom, start, end = condition.split('_')[:3]
            chrom = int(chrom[3:])
            start = int(start)
            end = int(end)
            outfname = f'{ukb}/association/plots/{phenotype}{run_suffix}.manhattan.{condition}.{ext}'
            conditioned_isnps = get_conditioned_isnps(condition)
            conditioned_strs = get_conditioned_strs(condition)

            unconditioned_str_results, unconditioned_snp_results, _, _ = load_data(
                phenotype, binary, None
            )

        if bool(args.peaks_spacing):
            print("Adding peak data ... ", end='', flush=True)
            start_time = time.time()
            # change the output file
            outfname = '.'.join(outfname.split('.')[:-1]) + f'.peaks_{args.peaks_spacing}_{args.peaks_thresh}.{ext}'

            peaks_fname = f'{ukb}/signals/peaks/{phenotype}_{args.peaks_spacing}_{args.peaks_thresh}.tab'
            peaks = pd.read_csv(
                peaks_fname,
                delimiter='\t',
                header=0,
                dtype=utils.get_dtypes(peaks_fname, {'ref_(snp_only)': object, 'alt_(snp_only)': object})
            )
            peaks.rename(columns={'chrom': 'chr', 'ref_(snp_only)': 'ref', 'alt_(snp_only)': 'alt'}, inplace=True)
            peaks['tagged_by_other_variant_type'] = peaks['tagged_by_other_variant_type'].astype(float)

            str_peaks = peaks[peaks['variant_type'] == 'STR']
            my_str_results = pd.DataFrame.from_records(my_str_results)
            my_str_results = my_str_results.merge(
                str_peaks[['chr', 'pos', 'tagged_by_other_variant_type']],
                how='left',
                on=['chr', 'pos']
            )
            my_str_results['is_peak'] = ~np.isnan(my_str_results['tagged_by_other_variant_type'])
            my_str_results = utils.df_to_recarray(my_str_results)

            snp_peaks = peaks[peaks['variant_type'] == 'SNP']
            plink_snp_results = pd.DataFrame.from_records(plink_snp_results)
            plink_snp_results = plink_snp_results.merge(
                snp_peaks[['chr', 'pos', 'ref', 'alt', 'tagged_by_other_variant_type']],
                how='left',
                on=['chr', 'pos', 'ref', 'alt']
            )
            plink_snp_results['is_peak'] = ~np.isnan(plink_snp_results['tagged_by_other_variant_type'])
            plink_snp_results = utils.df_to_recarray(plink_snp_results)
            print(f"done ({time.time() - start_time:.2e}s)", flush=True)
        elif args.finemap_signals:
            # change the output file
            outfname = '.'.join(outfname.split('.')[:-1]) + f'.FINEMAP.{ext}'
            finemap_loc = f'{ukb}/finemapping/finemap_results/{phenotype}'
            out_files = [f'{finemap_loc}/{d}/finemap_output.snp' for d in
                         os.listdir(finemap_loc) if d[0] in '0123456789']
            out_files = [f for f in out_files if os.path.exists(f) and os.path.getsize(f) > 0]
            snp_finemap_signals, str_finemap_signals, finemap_regions = load_finemap_signals(sorted(
                out_files,
                key=lambda dirname: [int(val) for val in dirname.split('/')[-2].split('_')]
            ))
    else:
        chrom = args.chrom
        start = args.start
        end = args.end
        gwas_catalog = gwas_catalog_ids = None
        my_str_results = my_str_run_date = None

        if not binary:
            out_runtype = 'continuous'
            readme_suffix = ''
        else:
            out_runtype = 'binary_' + binary
            readme_suffix = '_' + binary

        if not args.only_plink:
            my_snp_results = load_my_snp_results(phenotype, binary, chrom, start, end)
        else:
            outfname = f'{ukb}/association/plots/{out_runtype}_plink_snp.{ext}'
            my_snp_results = None

        if not args.only_mine:
            plink_snp_results = load_plink_results(phenotype, binary, None, fname=args.plink_snp_fname)
            plink_snp_results = plink_snp_results[plink_snp_results['chr'] == chrom]
        else:
            outfname = f'{ukb}/association/plots/{out_runtype}_my_imputed_snp.{ext}'
            plink_snp_results = None

        if not args.only_mine and not args.only_plink:
            outfname = f'{ukb}/association/plots/{out_runtype}_my_imputed_snp_vs_plink.{ext}'

        with open(f'{ukb}/association/results/{phenotype}/my_imputed_snp{readme_suffix}/README.txt') as README:
            date_line = next(README)
            my_snp_run_date = date_line.split(' ')[2]

    if not condition:
        # TODO for future runs uncomment
        '''
        with open(f'{ukb}/association/results/{phenotype}/plink_snp/logs/chr21.plink.stdout') as outlog:
            plink_snp_run_date = ' '.join(outlog.readlines()[-1].split(' ')[2:])
        '''
        plink_snp_run_date = None
    else:
        plink_snp_run_date = None

    with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt') as unit_file:
        unit = next(unit_file).strip()

    if not args.overview and not condition:
        make_manhattan_plots(
            outfname,
            phenotype,
            binary,
            unit,
            my_str_results,
            my_str_run_date,
            plink_snp_results,
            plink_snp_run_date,
            gwas_catalog,
            gwas_catalog_ids,
            my_snp_results,
            my_snp_run_date,
            chrom = chrom,
            start = start,
            end = end,
            snp_finemap_signals = snp_finemap_signals,
            str_finemap_signals = str_finemap_signals,
            finemap_regions = finemap_regions,
            conditioned_isnps = conditioned_isnps,
            conditioned_strs = conditioned_strs
        )
    elif args.overview:
        make_overview_manhattan(
            f'{ukb}/association/plots/{phenotype}.overview.manhattan.png',
            phenotype,
            my_str_results,
            plink_snp_results
        )
    else:
        conditioned_man = make_manhattan_plots(
            '.png',
            phenotype,
            binary,
            unit,
            my_str_results,
            my_str_run_date,
            plink_snp_results,
            plink_snp_run_date,
            gwas_catalog,
            gwas_catalog_ids,
            my_snp_results,
            my_snp_run_date,
            chrom = chrom,
            start = start,
            end = end,
            snp_finemap_signals = snp_finemap_signals,
            str_finemap_signals = str_finemap_signals,
            finemap_regions = finemap_regions,
            conditioned_isnps = conditioned_isnps,
            conditioned_strs = conditioned_strs,
            return_figure = True
        )
        unconditioned_man = make_manhattan_plots(
            '.png',
            phenotype,
            binary,
            unit,
            unconditioned_str_results,
            None,
            unconditioned_snp_results,
            None,
            gwas_catalog,
            gwas_catalog_ids,
            my_snp_results,
            my_snp_run_date,
            chrom = chrom,
            start = start,
            end = end,
            snp_finemap_signals = snp_finemap_signals,
            str_finemap_signals = str_finemap_signals,
            finemap_regions = finemap_regions,
            conditioned_isnps = None,
            conditioned_strs = None,
            return_figure = True
        )
        layout = bokeh.layouts.grid([unconditioned_man, conditioned_man])
        conditioned_man.y_range = unconditioned_man.y_range
        conditioned_man.title.text = 'Conditional ' + conditioned_man.title.text
        for man in (unconditioned_man, conditioned_man):
            man.background_fill_color = None
            man.border_fill_color = None
        bokeh.io.export_png(layout, filename=outfname)

if __name__ == "__main__":
    main()

