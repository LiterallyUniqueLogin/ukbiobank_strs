#!/usr/bin/env python3

import argparse
import copy
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
import region_plots

# python stops being able to distinguish between small numbers and zero at 1e-324 == 0
# so cutoff a bit before then
max_p_val = 300 # in -log10

def create_source_dict(
        data: np.recarray,
        start_chrom: Optional[int],
        cols_to_skip: Set[str] = set(),
        cols_to_include: Set[str] = set(),
        chrs_to_var_signals = None,
        start_display_cap = None,
        use_tstat = False
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
    if use_tstat:
        if 't_stat' not in data.dtype.names:
            print('t_stat', flush=True)
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
        if start_display_cap is not None:
            if not use_tstat:
                sources[chrom].data['display_val'] = np.minimum(
                    sources[chrom].data['p_val'], start_display_cap
                )
            else:
                sources[chrom].data['display_val'] = np.minimum(
                    sources[chrom].data['t_stat'], start_display_cap
                )
        else:
            if not use_tstat:
                sources[chrom].data['display_val'] = sources[chrom].data['p_val'].copy()
            else:
                sources[chrom].data['display_val'] = sources[chrom].data['t_stat'].copy()

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
        plink_snp_data,
        #gwas_catalog,
        #gwas_catalog_ids,
        my_snp_data,
        chrom,
        start,
        end,
        snp_finemap_signals,
        str_finemap_signals,
        finemap_regions,
        conditioned_strs,
        conditioned_isnps,
        *,
        return_figure = False,
        legend=True,
        chr_lens,
        use_tstat):
    ext = outfname.split('.')[-1]

    if not binary:
        stat = 'mean'
    else:
        stat = 'fraction'

    conditioning = conditioned_strs or conditioned_isnps

    if return_figure:
        assert outfname[0] == '.'

    plot_my_str_data = my_str_data is not None
    plot_my_snp_data = my_snp_data is not None
    plot_plink_snp_data = plink_snp_data is not None
    #plot_gwas_catalog = gwas_catalog is not None
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
        'total_best_guess_genotypes',
        'subset_total_best_guess_genotypes'
    }

    sources = []
    source_dicts = []

    if not conditioning and not start:
        if not use_tstat:
            start_display_cap = 30
        else:
            start_display_cap = 10 # TODO
    else:
        start_display_cap = None

    if plot_my_str_data:
        my_str_source, my_str_sources = create_source_dict(
            my_str_data,
            chrom,
            cols_to_skip = cols_to_skip,
            chrs_to_var_signals = str_finemap_signals,
            start_display_cap = start_display_cap,
            use_tstat = use_tstat,
        )
        sources.append(my_str_source)
        source_dicts.append(my_str_sources)
        for source in (my_str_source, *my_str_sources.values()):
            source.data['size'] = np.full(source.data['pos'].shape, 12)
            if plot_peaks:
                source.data['size'][source.data['is_peak']] = 20

    if plot_my_snp_data:
        my_snp_source, my_snp_sources = create_source_dict(
            my_snp_data,
            chrom,
            cols_to_skip = cols_to_skip,
            start_display_cap = start_display_cap,
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
            chrs_to_var_signals = snp_finemap_signals,
            use_tstat = use_tstat
        )
        sources.append(plink_snp_source)
        source_dicts.append(plink_snp_sources)
        for source in (plink_snp_source, *plink_snp_sources.values()):
            source.data['size'] = np.full(source.data['pos'].shape, 10)
            if plot_peaks:
                source.data['size'][source.data['is_peak']] = 20

#    if plot_gwas_catalog:
#        catalog_source, catalog_sources = create_source_dict(np.rec.fromarrays((
#            gwas_catalog[:, 0],
#            gwas_catalog[:, 1],
#            gwas_catalog[:, 2],
#            gwas_catalog_ids
#        ), names=['chr', 'pos', 'p_val', 'rsids']),
#        chrom)
#        sources.append(catalog_source)
#        source_dicts.append(catalog_sources)

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
        locus_plot.title.text_font_size = '18px'
        locus_plot.axis.axis_label_text_font_size = '18px'
        locus_plot.axis.major_label_text_font_size = '14px'

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

        locus_plot.varea('x', 'CI5e_2_upper', 'CI5e_8_upper', source=locus_source, color="red", alpha=0.2, legend_label='1 - 5e-8 Confidence')
        locus_plot.varea('x', 'CI5e_2_upper', 'CI5e_2_lower', source=locus_source, color="red", alpha=0.4, legend_label='0.95 confidence')
        locus_plot.varea('x', 'CI5e_2_lower', 'CI5e_8_lower', source=locus_source, color="red", alpha=0.2)
        locus_plot.line('x', 'y', source=locus_source, line_width=2, color="black")
        locus_plot.circle('x', 'y', source=locus_source, color="black", size=6, legend_label='mean')
        locus_plot.legend.label_text_font_size = '9px'
        
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
        x_axis_label = f'Mbp chr{chrom} (hg19)'

    if ext == 'html':
        output_backend = 'webgl'
    elif ext == 'svg':
        output_backend = 'svg'
    else:
        assert ext == 'png'
        output_backend = 'canvas'

    if start_display_cap is not None:
        start_height_cap = start_display_cap
    else:
        if plot_plink_snp_data:
            plink_data = plink_snp_source.data
            start_height_cap = np.max(
                plink_data['display_val'][(start <= plink_data['pos']) & (plink_data['pos'] <= end)]
            )

        if plot_my_str_data:
            str_data = my_str_source.data
            included_strs = str_data['display_val'][(start <= str_data['pos']) & (str_data['pos'] <= end)]
            if included_strs.shape[0] > 0:
                start_height_cap = max(start_height_cap, np.max(included_strs))

    manhattan_plot = bokeh.plotting.figure(
        width=1600,
        height=400,
        x_axis_label=x_axis_label,
        y_axis_label='-log10(p-value)',
        tools='xzoom_in,xzoom_out,save',
        y_range=(-0.025*start_height_cap, start_height_cap*1.025),
        output_backend=output_backend
    )
    manhattan_plot.axis.axis_label_text_font_size = '26px'
    manhattan_plot.axis.major_label_text_font_size = '20px'
    manhattan_plot.xaxis.formatter = bokeh.models.FuncTickFormatter(code='''
        return Math.floor(tick/100000)/10
    ''')

    if ext != 'html':
        manhattan_plot.grid.grid_line_color = None
        manhattan_plot.background_fill_color = None
        manhattan_plot.border_fill_color = None
        manhattan_plot.toolbar_location = None

    if start:
        #x_width = end - start
        #manhattan_plot.x_range = bokeh.models.Range1d(start-0.025*x_width, end+0.025*x_width)
        manhattan_plot.x_range = bokeh.models.Range1d(start, end)

    if ext == 'html':
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
        legend_label='genome-wide significance threshold'
    )
    manhattan_plot.legend.label_text_font_size = '22px'

    if finemap_regions:
        cmap_field_name = 'FINEMAP_pcausal'
    else:
        cmap_field_name = 'pos' # arbitrary existant field

    if plot_plink_snp_data:
        plink_snp_color = '#00B8FF'
        plink_snp_cmap = bokeh.transform.linear_cmap(
            field_name = cmap_field_name,
            low = 0,
            high = 0.2,
            palette = colorcet.kr,
            nan_color = 'purple'
        )
        plink_snp_manhattan = manhattan_plot.circle(
            'pos',
            'display_val',
            source=plink_snp_source,
            legend_label='SNPs and indels',
            color=plink_snp_color,
            size='size',
            muted_alpha=0.1
        )
    if plot_my_str_data:
        my_str_color = '#FF520D'
        my_str_cmap = bokeh.transform.linear_cmap(
            field_name = cmap_field_name,
            low = 0,
            high = 0.2,
            palette = colorcet.kr,
            nan_color = 'purple'
        )
        my_str_manhattan = manhattan_plot.square_pin(
            'pos',
            'display_val',
            source=my_str_source,
            legend_label='STRs',
            color=my_str_color,
            size='size',
            muted_alpha=0.1
        )
    if plot_my_snp_data:
        my_snp_manhattan = manhattan_plot.circle(
            'pos',
            'display_val',
            source=my_snp_source,
            legend_label='SNPs my code',
            color=(204, 121, 167),
            size=5,
            muted_alpha=0.1
        )
#    if plot_gwas_catalog:
#        catalog_manhattan = manhattan_plot.square(
#            'pos',
#            'display_val',
#            source=catalog_source,
#            legend_label='GWAS Catalog Hits',
#            color=(213, 94, 0),
#            size=7,
#            muted_alpha=0.1
#        )
    condition_color='#CD00E0'
    field_name = "p_val" if not use_tstat else "t_stat"
    if conditioned_isnps:
        manhattan_plot.circle(
            [snp[0] for snp in conditioned_isnps],
            [plink_snp_data[(plink_snp_data['pos'] == snp[0]) & (plink_snp_data['ref'] == snp[1]) & (plink_snp_data['alt'] == snp[2])][0][field_name] for snp in conditioned_isnps],
            legend_label='Conditioned-on variants',
            line_color='black',
            fill_color=None,
            size=20,
            line_width=3
        )
    if conditioned_strs:
        # window to account for flanking bps issue
        temp_datas = [my_str_data[(my_str_data['pos'] >= str_ - 10) & (my_str_data['pos'] <= str_ + 10)] for str_ in conditioned_strs]
        for temp_data, str_ in zip(temp_datas, conditioned_strs):
            assert temp_data.shape[0] == 1, (my_str_data['pos'], str_, temp_data)
        manhattan_plot.circle(
            conditioned_strs,
            [temp_data[0][field_name] for temp_data in temp_datas],
            legend_label='Conditioned-on variants',
            line_color='black',
            fill_color=None,
            size=20,
            line_width=3
        )

    if ext == 'html':
        # add hover tooltips for each manhattan plot
        # see https://stackoverflow.com/questions/49282078/multiple-hovertools-for-different-lines-bokeh
        hover_tools = []
        if plot_my_str_data:
            my_str_hover = bokeh.models.tools.HoverTool(renderers=[my_str_manhattan])
            my_str_hover.tooltips = [
                ('var type', 'STR'),
                ('alleles:', '@alleles'),
                ('pos', '@pos'),
                ('-log10(p_val) my code', '@p_val')
            ]
            if use_tstat:
                my_str_hover.tooltips.append(
                    ('t_stat', '@t_stat')
                )
            if binary:
                my_str_hover.tooltips.append(
                    ('Using Firth penalization during regression?', '@firth')
                )
            if conditioning:
                my_str_hover.tooltips.append(
                    ('-log10(unconditioned_p_val) my code', '@unconditioned_p')
                )
                if use_tstat:
                    my_str_hover.tooltips.append(
                        ('unconditioned_t_stat', '@unconditioned_t_stat')
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

        if plot_my_snp_data:
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

        if plot_plink_snp_data:
            plink_snp_hover = bokeh.models.tools.HoverTool(renderers=[plink_snp_manhattan])
            plink_snp_hover.tooltips = [
                ('var type', 'SNP'),
                ('alleles:', '@alleles'),
                ('pos', '@pos'),
                ('ID', '@id'),
                ('-log10(p_val) Plink', '@p_val')
            ]
            if use_tstat:
                my_str_hover.tooltips.append(
                    ('t_stat', '@t_stat')
                )
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
            if conditioning:
                plink_snp_hover.tooltips.append(
                    ('-log10(unconditioned_p_val) Plink', '@unconditioned_p')
                )
                if use_tstat:
                    my_str_hover.tooltips.append(
                        ('unconditioned_t_stat', '@unconditioned_t_stat')
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

#        if plot_gwas_catalog:
#            catalog_hover = bokeh.models.tools.HoverTool(renderers=[catalog_manhattan])
#            catalog_hover.tooltips = [
#                ('var type', 'GWAS Catalog SNP'),
#                ('pos', '@pos'),
#                ('rsid:', '@rsids'),
#                ('-log10(p_val)', '@p_val')
#            ]
#            manhattan_plot.add_tools(catalog_hover)
#            hover_tools.append(catalog_hover)

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

                locus_plot.title.text = 'STR ' + data['chr'][idx].toString() + ':' + data['pos'][idx].toString() + ' Repeat Unit: ' + data['motif'][idx];
            """
        )
        my_str_source.selected.js_on_change('indices', plot_tap_callback)

    height_slider = bokeh.models.Slider(
        start = 8,
        end=max_p_val if not use_tstat else 100, #TODO
        value=start_height_cap,
        step=1,
        title="p-value cap" if not use_tstat else "t_stat cap",
        sizing_mode = 'stretch_width'
    )
    for source in sources:
        height_callback = bokeh.models.CustomJS(
            args=dict(source=source, height_slider=height_slider, use_tstat = use_tstat),
            code="""
                const data = source.data;
                const display_val = data['display_val'];
                if (!use_tstat) {
                    const display_source= data['p_val'];
                } else {
                    const display_source = data['t_stat'];
                }
                const new_max = height_slider.value;
                for (var i = 0; i < display_val.length; i++) {
                    display_val[i] = Math.min(display_source[i], new_max);
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
                    range=manhattan_plot.x_range,
                    use_tstat = use_tstat,
                ),
                code = """
                    source.data = source_dict[this.value].data;
                    const data = source.data;
                    const display_val = data['display_val'];
                    if (!use_tstat) {
                        const display_source= data['p_val'];
                    } else {
                        const display_source = data['t_stat'];
                    }
                    const new_max = height_slider.value;
                    for (var i = 0; i < display_val.length; i++) {
                        display_val[i] = Math.min(display_source[i], new_max);
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

    manhattan_plot.legend.click_policy="mute"

    manhattan_plot.legend.visible = legend
    manhattan_plot.legend.location = 'top_left'
    '''
    if legend:
        manhattan_plot.add_layout(manhattan_plot.legend[0], 'right')
        del manhattan_plot.legend[0]
    '''
    if return_figure:
        return manhattan_plot

    if ext == 'html':
        html = bokeh.embed.file_html(layout, bokeh.resources.CDN, f'Manhattan plot {phenotype}')
        with open(outfname, 'w') as outfile:
            outfile.write(html)
    elif ext == 'png':
        manhattan_plot.background_fill_color = None
        manhattan_plot.border_fill_color = None
        bokeh.io.export_png(manhattan_plot, filename=outfname)
    else:
        assert ext == 'svg'
        bokeh.io.export_svg(manhattan_plot, filename=outfname)
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outfname')
    parser.add_argument('phenotype')
    parser.add_argument('unit')
    parser.add_argument('chr_lens_fname')
    parser.add_argument('str_gwas_results')
    parser.add_argument('snp_gwas_results')
    parser.add_argument('ext', help='file extension', choices=['png', 'html', 'svg'])
    #parser.add_argument("--my-plink-comparison", action='store_true', default=False)
    #parser.add_argument("--only-plink", action='store_true', default=False)
    #parser.add_argument("--only-mine", action='store_true', default=False)
    parser.add_argument("--chrom", type=int)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument('--conditioned-STRs', nargs='*', type=int)
    parser.add_argument('--conditioned-imputed-SNPs', nargs='*')
    parser.add_argument('--conditional-STR-results')
    parser.add_argument('--conditional-imputed-SNP-results')
    parser.add_argument('--finemap-output-files', nargs='+')
    parser.add_argument('--binary', default=False, choices={'linear', 'logistic'})
    #parser.add_argument('--plink-snp-fname')
    parser.add_argument('--peaks-fname')
    parser.add_argument('--use-tstat', action='store_true', default=False)
    args = parser.parse_args()

    phenotype = args.phenotype
    binary = args.binary
    chrom = args.chrom
    start = args.start
    end = args.end
    condition = bool(args.conditioned_imputed_SNPs) or bool(args.conditioned_STRs)
    chr_lens = np.genfromtxt(
        args.chr_lens_fname,
        skip_header=1,
        usecols=(1),
        dtype=int
    )

    assert condition == (args.conditional_STR_results is not None) == (args.conditional_imputed_SNP_results is not None)

    if not args.binary:
        run_suffix = ''
    else:
        run_suffix = '.' + binary

    # assert bool(args.finemap_output_files) + args.my_plink_comparison + condition + bool(args.peaks_fname) <= 1
    assert bool(args.finemap_output_files) + condition + bool(args.peaks_fname) <= 1

    assert bool(start) == bool(end)
    if bool(start):
        assert bool(chrom)

    if condition:
        assert bool(chrom) and bool(start)
        conditioned_isnps = []
        for SNP in args.conditioned_imputed_SNPs:
            pos, ref, alt = SNP.split('_')
            pos = int(pos)
            conditioned_isnps.append((pos, ref, alt))
        conditioned_strs = args.conditioned_STRs
    else:
        conditioned_isnps = conditioned_strs = None

    #if args.my_plink_comparison:
    #    assert bool(chrom)
    #if bool(args.plink_snp_fname) or args.only_mine or args.only_plink:
    #    assert args.my_plink_comparison
    #assert args.only_mine + args.only_plink <= 1

    if not bool(args.finemap_output_files):
        snp_finemap_signals = None
        str_finemap_signals = None
        finemap_regions = None

    if binary:
        runtype_suffix = '_' + binary
    else:
        runtype_suffix = ''

    #if not bool(args.my_plink_comparison):
    my_str_results = region_plots.load_my_str_results(
        phenotype, binary, args.str_gwas_results,
        args.conditional_STR_results
    )
    plink_snp_results = region_plots.load_plink_results(
        phenotype, binary, args.snp_gwas_results,
        args.conditional_imputed_SNP_results
    )
    #gwas_catalog, gwas_catalog_ids = region_plots.load_gwas_catalog(phenotype)

    my_snp_results = None
    if condition:
        unconditioned_str_results = region_plots.load_my_str_results(
            phenotype, binary, args.str_gwas_results
        )
        unconditioned_snp_results = region_plots.load_plink_results(
            phenotype, binary, args.snp_gwas_results
        )

    if bool(args.peaks_fname):
        print("Adding peak data ... ", end='', flush=True)
        start_time = time.time()

        #peaks_fname = f'{ukb}/signals/peaks/{phenotype}_{args.peaks_spacing}_{args.peaks_thresh}.tab'
        peaks = pd.read_csv(
            args.peaks_fname,
            delimiter='\t',
            header=0,
            dtype=utils.get_dtypes(args.peaks_fname, {'ref_(snp_only)': object, 'alt_(snp_only)': object})
        )
        peaks.rename(columns={'chrom': 'chr', 'ref_(snp_only)': 'ref', 'alt_(snp_only)': 'alt'}, inplace=True)
        peaks['marker'] = 1

        str_peaks = peaks[peaks['variant_type'] == 'STR']
        my_str_results = pd.DataFrame.from_records(my_str_results)
        my_str_results = my_str_results.merge(
            str_peaks[['chr', 'pos', 'marker']],
            how='left',
            on=['chr', 'pos']
        )
        my_str_results['is_peak'] = ~np.isnan(my_str_results['marker'])
        my_str_results = utils.df_to_recarray(my_str_results)

        snp_peaks = peaks[peaks['variant_type'] == 'SNP']
        plink_snp_results = pd.DataFrame.from_records(plink_snp_results)
        plink_snp_results = plink_snp_results.merge(
            snp_peaks[['chr', 'pos', 'ref', 'alt', 'marker']],
            how='left',
            on=['chr', 'pos', 'ref', 'alt']
        )
        plink_snp_results['is_peak'] = ~np.isnan(plink_snp_results['marker'])
        plink_snp_results = utils.df_to_recarray(plink_snp_results)
        print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    elif bool(args.finemap_output_files):
        """
        finemap_loc = f'{ukb}/finemapping/finemap_results/{phenotype}'
        out_files = []
        with open(f'{ukb}/signals/regions/{phenotype}.tab') as regions:
            next(regions)
            for line in regions:
                region_chrom, region_start, region_end, any_strs = line.strip().split('\t')
                if any_strs == 'False':
                    continue
                out_files.append(f'{finemap_loc}/{region_chrom}_{region_start}_{region_end}/finemap_output.snp')
        """
        snp_finemap_signals, str_finemap_signals, finemap_regions = region_plots.load_finemap_signals(args.finemap_output_files)
#    else: # my plink comparison
#        gwas_catalog = gwas_catalog_ids = None
#
#        if not binary:
#            out_runtype = 'continuous'
#            readme_suffix = ''
#        else:
#            out_runtype = 'binary_' + binary
#            readme_suffix = '_' + binary
#
#        if not args.only_plink:
#            my_snp_results = region_plots.load_my_snp_results(
#                phenotype, binary, chrom, start, end
#            )
#        else:
#            my_snp_results = None
#
#        if not args.only_mine:
#            plink_snp_results = region_plots.load_plink_results(
#                phenotype, binary, args.plink_snp_fname
#            )
#            plink_snp_results = plink_snp_results[plink_snp_results['chr'] == chrom]
#        else:
#            plink_snp_results = None

    if not condition:
        make_manhattan_plots(
            args.outfname,
            phenotype,
            binary,
            args.unit,
            my_str_results,
            plink_snp_results,
            #gwas_catalog,
            #gwas_catalog_ids,
            my_snp_results,
            chrom = chrom,
            start = start,
            end = end,
            snp_finemap_signals = snp_finemap_signals,
            str_finemap_signals = str_finemap_signals,
            finemap_regions = finemap_regions,
            conditioned_isnps = conditioned_isnps,
            conditioned_strs = conditioned_strs,
            chr_lens = chr_lens,
            use_tstat = args.use_tstat,
        )
    else:
        conditioned_man = make_manhattan_plots(
            f'.{args.ext}',
            phenotype,
            binary,
            args.unit,
            my_str_results,
            plink_snp_results,
            #gwas_catalog,
            #gwas_catalog_ids,
            my_snp_results,
            chrom = chrom,
            start = start,
            end = end,
            snp_finemap_signals = snp_finemap_signals,
            str_finemap_signals = str_finemap_signals,
            finemap_regions = finemap_regions,
            conditioned_isnps = None,
            conditioned_strs = None,
            return_figure = True,
            legend=False,
            chr_lens = chr_lens,
            use_tstat = args.use_tstat,
        )
        unconditioned_man = make_manhattan_plots(
            f'.{args.ext}',
            phenotype,
            binary,
            args.unit,
            unconditioned_str_results,
            unconditioned_snp_results,
            #gwas_catalog,
            #gwas_catalog_ids,
            my_snp_results,
            chrom = chrom,
            start = start,
            end = end,
            snp_finemap_signals = snp_finemap_signals,
            str_finemap_signals = str_finemap_signals,
            finemap_regions = finemap_regions,
            conditioned_isnps = conditioned_isnps,
            conditioned_strs = conditioned_strs,
            return_figure = True,
            chr_lens = chr_lens,
            use_tstat = args.use_tstat,
        )
        layout = bokeh.layouts.grid([unconditioned_man, conditioned_man])
        conditioned_man.x_range = unconditioned_man.x_range
        conditioned_man.y_range = unconditioned_man.y_range
        for man in (unconditioned_man, conditioned_man):
            man.background_fill_color = None
            man.border_fill_color = None

        if args.ext == 'html':
            html = bokeh.embed.file_html(layout, bokeh.resources.CDN, 'conditional manhattan')
            with open(args.outfname, 'w') as outfile:
                outfile.write(html)
        else:
            for man in unconditioned_man, conditioned_man:
                man.toolbar_location = None
            if args.ext == 'png':
                bokeh.io.export_png(layout, filename=args.outfname)
            else:
                assert args.ext == 'svg'
                bokeh.io.export_svg(layout, filename=args.outfname)

if __name__ == "__main__":
    main()

