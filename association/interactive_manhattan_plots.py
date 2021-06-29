#!/usr/bin/env python3

import argparse
import ast
import copy
import os
import os.path
import time
from typing import Dict, Tuple, Optional, Set

import bokeh.colors
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

cum_lens = np.cumsum(chr_lens)

start_p_val_cap = 30
max_p_val = 350

def create_source_dict(
        data: np.recarray,
        start_chrom: Optional[int],
        cols_to_skip: Set[str] = set(),
        cols_to_include: Set[str] = set(),
        chrs_to_var_signals = None
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
            1/0
    for chrom in cds_range:
        #print(f'Converting chrom {chrom} ...', end='\r')
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
        sources[chrom].data['display_p_val'] = np.minimum(
            sources[chrom].data['p_val'], start_p_val_cap
        )
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
            sources[chrom].data['FINEMAP_pcausal'] = merged_data['p']

    copy_source = bokeh.models.ColumnDataSource(copy.deepcopy(sources[start_chrom].data))

    return copy_source, sources

def make_manhattan_plots(
        outfname,
        phenotype,
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
        finemap_regions):

    plot_my_str_data = my_str_data is not None
    plot_my_snp_data = my_snp_data is not None
    plot_gwas_catalog = gwas_catalog is not None

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

    if plot_my_str_data:
        my_str_source, my_str_sources = create_source_dict(
            my_str_data,
            chrom,
            cols_to_skip = cols_to_skip,
            chrs_to_var_signals = str_finemap_signals
        )
        sources.append(my_str_source)
        source_dicts.append(my_str_sources)

    if plot_my_snp_data:
        my_snp_source, my_snp_sources = create_source_dict(
            my_snp_data,
            chrom,
            cols_to_skip = cols_to_skip,
        )
        sources.append(my_snp_source)
        source_dicts.append(my_snp_sources)

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

    if plot_my_str_data:
        locus_plot = bokeh.plotting.figure(
            width=400,
            height=400,
            title='Click on an STR in the Manhattan plot!                                      ',
            x_axis_label='Sum of allele lengths (repeat copies)',
            y_axis_label=f'Mean {phenotype} ({unit})',
            tools='save',
            y_range=(-2, 5)
        )
        locus_plot.add_layout(bokeh.models.Title(
            text="Phenotype values are unadjusted for covariates or conditioned genotypes",
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
        x_axis_label = 'Position'
    else:
        x_axis_label = f'Position (chr {chrom})'

    manhattan_plot = bokeh.plotting.figure(
        width=1200,
        height=900,
        title=(phenotype.capitalize() + ' Manhattan Plot'),
        x_axis_label=x_axis_label,
        y_axis_label='-log10(p-value)',
        tools='xzoom_in,xzoom_out,save',
        y_range=(-.25, start_p_val_cap+.25),
        output_backend="webgl"
    )

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

    plink_snp_color = bokeh.colors.RGB(0, 114, 178)
    plink_snp_cmap = bokeh.transform.linear_cmap(
        field_name = cmap_field_name,
        low = 0,
        high = 0.2,
        palette = colorcet.kr,
        nan_color = 'purple'
    )
    ''' 
        palette = [plink_snp_color],
        nan_color = plink_snp_color,
        low_color = plink_snp_color,
        high_color = plink_snp_color
    )
    '''
    plink_snp_manhattan = manhattan_plot.circle(
        'pos',
        'display_p_val',
        source=plink_snp_source,
        legend_label='SNPs Plink',
        color=plink_snp_color,
        size=5,
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
        '''
            palette = [my_str_color],
            nan_color = my_str_color,
            low_color = my_str_color,
            high_color = my_str_color
        )
        '''
        my_str_manhattan = manhattan_plot.square_pin(
            'pos',
            'display_p_val',
            source=my_str_source,
            legend_label='STRs my code',
            color=my_str_color,
            size=6
        )
    if plot_my_snp_data:
        my_snp_manhattan = manhattan_plot.circle(
            'pos',
            'display_p_val',
            source=my_snp_source,
            legend_label='SNPs my code',
            color=(204, 121, 167),
            size=5
        )
    if plot_gwas_catalog:
        catalog_manhattan = manhattan_plot.square(
            'pos',
            'display_p_val',
            source=catalog_source,
            legend_label='GWAS Catalog Hits',
            color=(213, 94, 0),
            size=7
        )

    # add custom tick formatter. See FuncTickFormatter
    # here: https://docs.bokeh.org/en/latest/docs/user_guide/styling.html#userguide-styling-axes-tick-label-formats
    setup_cum_len_js = "var cum_lens = [];"
    for cum_len in cum_lens:
        setup_cum_len_js += f"cum_lens.push({cum_len - chr_lens[0]});"

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
        if str_finemap_signals:
            my_str_hover.tooltips.append(('FINEMAP_pcausal', '@FINEMAP_pcausal{safe}'))
        str_means_start_idx = list(my_str_data.dtype.names).index(
            f'mean_{phenotype}_per_single_dosage'
        )
        for detail_name in my_str_data.dtype.names[7:str_means_start_idx]:
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
        snp_means_start_idx = list(my_snp_data.dtype.names).index(
            f'mean_{phenotype}_per_single_dosage'
        )
        for detail_name in my_snp_data.dtype.names[7:snp_means_start_idx]:
            if detail_name in cols_to_skip:
                continue
            my_snp_hover.tooltips.append((detail_name, f'@{detail_name}' '{safe}'))
        manhattan_plot.add_tools(my_snp_hover)
        hover_tools.append(my_snp_hover)

    plink_snp_hover = bokeh.models.tools.HoverTool(renderers=[plink_snp_manhattan])
    plink_snp_hover.tooltips = [
        ('var type', 'SNP'),
        ('alleles:', '@alleles'),
        ('pos', '@pos'),
        ('ID', '@id'),
        ('-log10(p_val) Plink', '@p_val'),
        ('Minor allele frequency', '@maf'),
        ('Imputation INFO', '@info')
    ]
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

    if plot_my_str_data:
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
                var locus_dict = JSON.parse(data['mean_{phenotype}_per_single_dosage'][idx]);
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
        value=start_p_val_cap,
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
            y_range.start = -.25
            y_range.end = new_max + 0.25;
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
            text=f"New plink SNP run dat: {plink_snp_run_date}",
            align='right'
        ), 'below')

    manhattan_plot.legend.click_policy="mute"

    html = bokeh.embed.file_html(layout, bokeh.resources.CDN, f'Manhattan plot {phenotype}')
    with open(outfname, 'w') as outfile:
        outfile.write(html)
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

def load_data(phenotype, condition):
    if not condition:
        return (
            load_my_str_results(phenotype, None),
            load_plink_results(phenotype, None),
            *load_gwas_catalog(phenotype)
        )
    else:
        return (
            load_my_str_results(phenotype, condition),
            load_plink_results(phenotype, condition),
            None,
            None
        )

def load_my_str_results(phenotype, condition):
    print(f"Loading my STR results for {phenotype} ... ", end='', flush=True)
    start_time = time.time()
    if not condition:
        str_results_fname = \
            f'{ukb}/association/plots/input/{phenotype}/my_str_results.tab'
    else:
        str_results_fname = \
            f'{ukb}/association/plots/input/{phenotype}/my_str_conditional_{condition}_results.tab'
    my_str_results = utils.df_to_recarray(pd.read_csv(
        str_results_fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        dtype=utils.get_dtypes(str_results_fname, {'locus_filtered': str})
    ))
    names = list(my_str_results.dtype.names)
    for idx, name in my_results_rename.items():
        names[idx] = name
    for idx, name in my_str_results_rename.items():
        names[idx] = name
    for colname in ('total_per_allele_dosages', 'total_hardcall_alleles',
                'subset_total_per_allele_dosages', 'subset_total_hardcall_alleles',
                'subset_allele_dosage_r2'):
        # convert allele lens from strings to floats, in addition round allele lens and values, but not NaN values
        new_col = np.array(list(map(
            lambda dict_str: {round(float(allele_len), 2): (round(val, 2) if val != 'NaN' else val) for allele_len, val in ast.literal_eval(dict_str).items()},
            my_str_results[colname]
        )))
        # convert allele_lens to ints if they are close enough
        new_col = np.array(list(map(
            lambda d: str({(int(key) if key == int(key) else key) : val for key, val in d.items()}),
            new_col
        )))
        my_str_results[colname] = new_col
    my_str_results.dtype.names = names
    my_str_results['p_val'] = -np.log10(my_str_results['p_val'])
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return my_str_results

def load_my_chr21_height_snp_results():
    print("Loading my chr21 height results ... ", end='', flush=True)
    start_time = time.time()
    snp_results_fname = f'{ukb}/association/plots/input/height/my_imputed_snp_chr21_results.tab'
    my_snp_results = utils.df_to_recarray(pd.read_csv(
        snp_results_fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        dtype=utils.get_dtypes(snp_results_fname, {'locus_filtered': str})
    ))
    names = list(my_snp_results.dtype.names)
    for idx, name in my_results_rename.items():
        names[idx] = name
    my_snp_results.dtype.names = names
    my_snp_results['p_val'] = -np.log10(my_snp_results['p_val'])
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return my_snp_results

def load_plink_results(phenotype, condition):
    # Load plink SNP results
    print(f"Loading plink SNP results for {phenotype} ... ", end='', flush=True)
    start_time = time.time()
    if not condition:
        plink_snp_fname = \
            f'{ukb}/association/plots/input/{phenotype}/plink_snp_results_with_mfi.npy'
    else:
        plink_snp_fname = \
            f'association/plots/input/{phenotype}/plink_snp_conditional_{condition}_with_mfi.npy'

    plink_results = np.load(plink_snp_fname)
    plink_results['p_val'] = -np.log10(plink_results['p_val'])

    assert np.all(plink_results['error'] == '.')

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return plink_results

def load_gwas_catalog(phenotype):
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
        # TODO uncomment once I've fixed the overlapping interavls bug
        '''
        chrs_to_snp_signals[i] = []
        chrs_to_str_signals[i] = []
        '''
        chrs_to_snp_signals[i] = {}
        chrs_to_str_signals[i] = {}
        chrs_to_regions[i] = []
    for signal in finemap_signals:
        region = signal.split('/')[-2]
        chrom, start, end = (int(val) for val in region.split('_'))
        chrs_to_regions[chrom].append(f'{start}_{end}')
        with open(signal) as per_var_output:
            next(per_var_output)
            for line in per_var_output:
                _id, pos, p = np.array(line.split())[[1, 3, 10]]
                pos = int(pos)
                split = _id.split('_')
                assert int(split[1]) == pos
                if p == 'NA':
                    p = np.nan
                else:
                    p = float(p)
                if _id[:4] == 'STR_':
                    '''
                    chrs_to_str_signals[chrom].append((pos, p))
                    '''
                    if pos not in chrs_to_str_signals[chrom] or chrs_to_str_signals[chrom][pos] > p:
                        chrs_to_str_signals[chrom][pos] = p
                elif _id[:4] == 'SNP_':
                    '''
                    chrs_to_snp_signals[chrom].append(
                        (pos, split[2], split[3], p)
                    )
                    '''
                    key = (pos, split[2], split[3])
                    if key not in chrs_to_snp_signals[chrom] or chrs_to_snp_signals[chrom][key] > p:
                        chrs_to_snp_signals[chrom][key] = p
                else:
                    raise ValueError(f'Found uninterpretable id {_id}')

    '''
    chrs_to_snp_signals = {
        key: pd.DataFrame(np.stack(val), columns=('pos', 'ref', 'alt', 'p'))
        for key, val in chrs_to_snp_signals.items()
    }
    chrs_to_str_signals = {
        key: pd.DataFrame(np.stack(val), columns=('pos', 'p'))
        for key, val in chrs_to_str_signals.items()
    }
    '''
    chrs_to_snp_signals = {
        key: pd.DataFrame.from_dict(val, orient='index', columns=['p'])
        for key, val in chrs_to_snp_signals.items()
    }
    for key, df in chrs_to_snp_signals.items():
        df.reset_index(inplace=True)
        df[['pos', 'ref', 'alt']] = pd.DataFrame(df['index'].tolist(), index=df.index)
        chrs_to_snp_signals[key] = df[['pos', 'ref', 'alt', 'p']]

    chrs_to_str_signals = {
        key: pd.DataFrame.from_dict(val, orient='index', columns=['p'])
        for key, val in chrs_to_str_signals.items()
    }
    for key, df in chrs_to_str_signals.items():
        df.reset_index(inplace=True)
    chrs_to_str_signals = {
        key: df.rename(columns={'index': 'pos'}) for key, df in chrs_to_str_signals.items()
    }

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return (chrs_to_snp_signals, chrs_to_str_signals, chrs_to_regions)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--my-plink-comparison", action='store_true')
    parser.add_argument("--phenotype")
    parser.add_argument("--condition")
    parser.add_argument('--finemap-signals', action='store_true', default=False)
    args = parser.parse_args()
    assert args.finemap_signals + bool(args.my_plink_comparison) + bool(args.condition) <= 1
    if not bool(args.my_plink_comparison):
        assert bool(args.phenotype)
    else:
        assert bool(args.phenotype)

    if not args.finemap_signals:
        snp_finemap_signals = None
        str_finemap_signals = None
        finemap_regions = None

    if not bool(args.my_plink_comparison):
        phenotype = args.phenotype
        my_str_results, plink_snp_results, gwas_catalog, gwas_catalog_ids = load_data(
            phenotype, args.condition
        )
        my_snp_results = None
        my_snp_run_date = None
        if not args.condition:
            chrom = None
            start = None
            end = None
            with open(f'{ukb}/association/results/{phenotype}/my_str/README.txt') as README:
                date_line = next(README)
                my_str_run_date = date_line.split(' ')[2]
            outfname = f'{ukb}/association/plots/{phenotype}_interactive_manhattan.html'
        else:
            my_str_run_date = None
            chrom, start, end = args.condition.split('_')[:3]
            chrom = int(chrom[3:])
            start = int(start)
            end = int(end)
            outfname = f'{ukb}/association/plots/{phenotype}_interactive_manhattan_c_{args.condition}.html'

        if args.finemap_signals:
            # change the output file
            outfname = f'{ukb}/association/plots/{phenotype}_interactive_manhattan_FINEMAP.html'
            finemap_loc = f'{ukb}/finemapping/finemap_results/overlapping_{phenotype}'
            #finemap_loc = f'{ukb}/finemapping/finemap_results/{phenotype}'
            out_files = [f'{finemap_loc}/{d}/finemap_output.snp' for d in
                         os.listdir(finemap_loc) if d[0] in '0123456789']
            out_files = [f for f in out_files if os.path.exists(f) and os.path.getsize(f) > 0]
            snp_finemap_signals, str_finemap_signals, finemap_regions = load_finemap_signals(sorted(
                out_files,
                key=lambda dirname: [int(val) for val in dirname.split('/')[-2].split('_')]
            ))
    else:
        phenotype = 'height'
        chrom = 21
        start = None
        end = None
        my_snp_results = load_my_chr21_height_snp_results()
        plink_snp_results = load_plink_results(phenotype, None)
        plink_snp_results = plink_snp_results[plink_snp_results['chr'] == 21]
        gwas_catalog = gwas_catalog_ids = None
        my_str_results = my_str_run_date = None

        outfname = f'{ukb}/association/plots/height_my_imputed_snp_vs_plink.html'
        with open(f'{ukb}/association/results/{phenotype}/my_str/README.txt') as README:
            date_line = next(README)
            my_snp_run_date = date_line.split(' ')[2]

    if not args.condition:
        with open(f'{ukb}/association/results/{phenotype}/plink_snp/logs/chr21.plink.stdout') as README:
            plink_snp_run_date = ' '.join(README.readlines()[-1].split(' ')[2:])
    else:
        plink_snp_run_date = None

    with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt') as unit_file:
        unit = next(unit_file).strip()

    make_manhattan_plots(
        outfname,
        phenotype,
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
        finemap_regions = finemap_regions
    )

if __name__ == "__main__":
    main()

