#!/usr/bin/env python3

import argparse
import copy
import os
import os.path
import time
from typing import Dict, Tuple, Optional, Set

import bokeh.embed
import bokeh.events
import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.models.callbacks
import bokeh.models.tools
import bokeh.plotting
import bokeh.resources
import numpy as np
import numpy.lib.recfunctions
import numpy.ma
import numpy.random
import pandas as pd

ukb = os.environ['UKB']

chr_lens = np.genfromtxt(
    f'{ukb}/misc_data/genome/chr_lens.txt',
    skip_header=1,
    usecols=(1),
    dtype=int
)

cum_lens = np.cumsum(chr_lens)

# from https://stackoverflow.com/questions/52579601/convert-dataframe-to-a-rec-array-and-objects-to-strings
def df_to_recarray(df):
    names = df.columns
    arrays = [df[col].values for col in names]

    formats = [ array.dtype if array.dtype != 'O'
                else f'{array.astype(str).dtype}' for array in arrays ]

    rec_array = np.rec.fromarrays(
        arrays,
        dtype={'names': names, 'formats': formats}
    )

    return rec_array

start_p_val_cap = 30
max_p_val = 350

def create_source_dict(
        data: np.recarray,
        cols_to_skip: Set[str] = set(),
        cols_to_include: Set[str] = set()
    ) -> Tuple[bokeh.models.ColumnDataSource, Dict[int, bokeh.models.ColumnDataSource]]:
    sources = {}
    assert len(cols_to_skip) == 0 or len(cols_to_include) == 0
    for field in 'chr', 'pos', 'p_val':
        if field not in data.dtype.names:
            print(field, flush=True)
            1/0
    for chrom in range(1, 23):
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
    copy_source_1 = bokeh.models.ColumnDataSource(copy.deepcopy(sources[1].data))

    return copy_source_1, sources

def plot_manhattan(plot, source, label, color, size=4):
    return plot.circle(
        'pos',
        'display_p_val',
        legend_label=label,
        source=source,
        color=color,
        size=size,
        muted_alpha=0.1
    )


def make_manhattan_plots(
        phenotype,
        unit,
        my_str_data,
        my_str_run_date,
        plink_snp_data,
        plink_snp_run_date,
        gwas_catalog,
        gwas_catalog_ids,
        my_snp_data,
        my_snp_run_date):

    plot_my_str_data = my_str_data is not None
    plot_my_snp_data = my_snp_data is not None
    plot_gwas_catalog = gwas_catalog is not None

    print(f"Plotting phenotype {phenotype} ... ", end='', flush=True)
    start_time = time.time()
    # reformat data for plotting

    cols_to_skip = {
        'coeff_phenotype',
        'coeff_intercept',
        'R^2',
        'total_hardcall_genotypes',
        'subset_total_hardcall_genotypes',
    }

    sources = []
    source_dicts = []

    if plot_my_str_data:
        my_str_source, my_str_sources = create_source_dict(my_str_data, cols_to_skip=cols_to_skip)
        sources.append(my_str_source)
        source_dicts.append(my_str_sources)

    if plot_my_snp_data:
        my_snp_source, my_snp_sources = create_source_dict(my_snp_data, cols_to_skip=cols_to_skip)
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
    plink_snp_source, plink_snp_sources = create_source_dict(plink_snp_data)
    sources.append(plink_snp_source)
    source_dicts.append(plink_snp_sources)

    if plot_gwas_catalog:
        catalog_source, catalog_sources = create_source_dict(np.rec.fromarrays((
            gwas_catalog[phenotype][:, 0],
            gwas_catalog[phenotype][:, 1],
            gwas_catalog[phenotype][:, 2],
            gwas_catalog_ids[phenotype]
        ), names=['chr', 'pos', 'p_val', 'rsids']))
        sources.append(catalog_source)
        source_dicts.append(catalog_sources)

    if plot_my_str_data:
        locus_plot = bokeh.plotting.figure(
            width=400,
            height=400,
            title='Click on an STR in the Manhattan plot!                                      ',
            x_axis_label='Sum of allele lengths (repeat copies)',
            y_axis_label=f'Mean {phenotype} ({unit})',
            tools='save'
        )

        locus_source = bokeh.models.ColumnDataSource({
            'x': [1, 2],
            'y': [1, 2],
            'CI5e_2': ['nan', 'nan'],
            'CI5e_8': ['nan', 'nan']
        })

        phenotype_means = locus_plot.circle('x', 'y', source=locus_source)

        means_hover = bokeh.models.tools.HoverTool(renderers=[phenotype_means])
        means_hover.tooltips = [
            ('5e-2 CI', '@CInorm'),
            ('5e-8 CI', '@CIgwas')
        ]
        locus_plot.add_tools(means_hover)

    # set up drawing canvas
    manhattan_plot = bokeh.plotting.figure(
        width=1200,
        height=900,
        title=(phenotype.capitalize() + ' Manhattan Plot'),
        x_axis_label='Position',
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

    line_source = bokeh.models.ColumnDataSource(dict(
        x=[0, chr_lens[1]],
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
    plink_snp_manhattan = plot_manhattan(
        manhattan_plot, plink_snp_source, 'SNPs Plink', color=(0, 114, 178)
    )
    if plot_my_str_data:
        my_str_manhattan = plot_manhattan(
            manhattan_plot, my_str_source, 'STRs my code', color=(204, 121, 167), size=6
        )
    if plot_my_snp_data:
        my_snp_manhattan = plot_manhattan(
            manhattan_plot, my_snp_source, 'SNPs my code', color=(204, 121, 167)
        )
    if plot_gwas_catalog:
        catalog_manhattan = plot_manhattan(
            manhattan_plot, catalog_source, 'GWAS Catalog Hits', color=(213, 94, 0), size=7
        )

    # add custom tick formatter. See FuncTickFormatter
    # here: https://docs.bokeh.org/en/latest/docs/user_guide/styling.html#userguide-styling-axes-tick-label-formats
    setup_cum_len_js = "var cum_lens = [];"
    for cum_len in cum_lens:
        setup_cum_len_js += f"cum_lens.push({cum_len - chr_lens[0]});"

    # add hover tooltips for each manhattan plot
    # see https://stackoverflow.com/questions/49282078/multiple-hovertools-for-different-lines-bokeh
    if plot_my_str_data:
        my_str_hover = bokeh.models.tools.HoverTool(renderers=[my_str_manhattan])
        my_str_hover.tooltips = [
            ('var type', 'STR'),
            ('alleles:', '@alleles'),
            ('pos', '@pos'),
            ('-log10(p_val) my code', '@p_val')
        ]
        str_means_start_idx = list(my_str_data.dtype.names).index(
            f'mean_residual_{phenotype}_per_single_dosage'
        )
        for detail_name in my_str_data.dtype.names[7:str_means_start_idx]:
            if detail_name in cols_to_skip:
                continue
            my_str_hover.tooltips.append((detail_name, f'@{detail_name}' '{safe}'))
        manhattan_plot.add_tools(my_str_hover)

    if plot_my_snp_data:
        my_snp_hover = bokeh.models.tools.HoverTool(renderers=[my_snp_manhattan])
        my_snp_hover.tooltips = [
            ('var type', 'STR'),
            ('alleles:', '@alleles'),
            ('pos', '@pos'),
            ('-log10(p_val) my code', '@p_val')
        ]
        snp_means_start_idx = list(my_snp_data.dtype.names).index(
            f'mean_residual_{phenotype}_per_single_dosage'
        )
        for detail_name in my_snp_data.dtype.names[7:snp_means_start_idx]:
            if detail_name in cols_to_skip:
                continue
            my_snp_hover.tooltips.append((detail_name, f'@{detail_name}' '{safe}'))
        manhattan_plot.add_tools(my_snp_hover)

    plink_snp_hover = bokeh.models.tools.HoverTool(renderers=[plink_snp_manhattan])
    plink_snp_hover.tooltips = [
        ('var type', 'SNP'),
        ('alleles:', '@alleles'),
        ('pos', '@pos'),
        ('ID', '@id'),
        ('-log10(p_val) Plink', '@p_val')
    ]
    manhattan_plot.add_tools(plink_snp_hover)

    if plot_gwas_catalog:
        catalog_hover = bokeh.models.tools.HoverTool(renderers=[catalog_manhattan])
        catalog_hover.tooltips = [
            ('var type', 'GWAS Catalog SNP'),
            ('pos', '@pos'),
            ('rsid:', '@rsids'),
            ('-log10(p_val)', '@p_val')
        ]
        manhattan_plot.add_tools(catalog_hover)

    if plot_my_str_data:
        manhattan_plot.toolbar.active_inspect = [my_str_hover, plink_snp_hover, catalog_hover]
    else:
        manhattan_plot.toolbar.active_inspect = [my_snp_hover, plink_snp_hover]

    if plot_my_str_data:
        tap = bokeh.models.tools.TapTool()
        manhattan_plot.add_tools(tap)
        manhattan_plot.toolbar.active_tap = tap

        plot_tap_callback = bokeh.models.CustomJS(
            args=dict(locus_plot=locus_plot, my_str_source=my_str_source, locus_source=locus_source),
            #args=dict(my_str_source=my_str_source, locus_source=locus_source),
            code = f"""
                if (my_str_source.selected.indices.length == 0) {{
                    return;
                }}
                var idx = my_str_source.selected.indices[0];
                var data = my_str_source.data;
                var locus_dict = JSON.parse(data['mean_residual_{phenotype}_per_single_dosage'][idx]);
                var CI5e_2_str = data['CI5e_2SingleDosagePhenotype'][idx].replaceAll('(', '[').replaceAll(')', ']').replaceAll('nan', '"NaN"');
                var CI5e_2_dict = JSON.parse(CI5e_2_str);
                var CI5e_8_str = data['CI5e_8SingleDosagePhenotype'][idx].replaceAll('(', '[').replaceAll(')', ']').replaceAll('nan', '"NaN"');
                var CI5e_8_dict = JSON.parse(CI5e_8_str);
                console.log(data['mean_residual_{phenotype}_per_single_dosage'][idx]);
                console.log(CI5e_2_str);
                console.log(CI5e_8_str);
                var gts = [];
                var str_gts = [];
                var phenotypes = [];
                var CI5e_2 = [];
                var CI5e_8 = []; 
                Object.entries(locus_dict).forEach(([key, value]) => {{
                    str_gts.push(key);
                    gts.push(parseFloat(key));
                    phenotypes.push(parseFloat(value));
                }});
                console.log(gts);
                for (idx=0; idx<gts.length; idx++) {{
                    console.log(idx, gts[idx], str_gts[idx]);
                    CI5e_2.push(CI5e_2_dict[str_gts[idx]].toString()); 
                    CI5e_8.push(CI5e_8_dict[str_gts[idx]].toString()); 
                }}
                
                var new_dict = {{'x': gts, 'y': phenotypes, 'CInorm': CI5e_2, 'CIgwas': CI5e_8}};
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

    if plot_my_str_data:
        layout = bokeh.layouts.grid([
            [None, locus_plot],
            [height_slider, chrom_select],
            manhattan_plot
        ])
    else:
        layout = bokeh.layouts.grid([
            [height_slider, chrom_select],
            manhattan_plot
        ])

    if plot_my_str_data:
        manhattan_plot.add_layout(bokeh.models.Title(
            text=f"My STR code run date: {my_str_run_date}",
            align='right'
        ), 'below')
    if plot_my_snp_data:
        manhattan_plot.add_layout(bokeh.models.Title(
            text=f"My SNP code run date: {my_snp_run_date}",
            align='right'
        ), 'below')
    manhattan_plot.add_layout(bokeh.models.Title(
        text=f"New plink SNP run dat: {plink_snp_run_date}",
        align='right'
    ), 'below')

    manhattan_plot.legend.click_policy="mute"

    html = bokeh.embed.file_html(layout, bokeh.resources.CDN, 'Manhattan plot test')
    if plot_my_snp_data:
        with open(f'{ukb}/association/plots/height_my_imputed_snp_vs_plink.html', 'w') as out_file:
            out_file.write(html)
    else:
        with open(f'{ukb}/association/plots/{phenotype}_interactive_manhattan.html', 'w') as out_file:
            out_file.write(html)
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

def get_dtypes(fname):
    df = pd.read_csv(
        fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        nrows=1
    )
    dtypes = dict(df.dtypes)
    dtypes['locus_filtered'] = str
    return dtypes

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

def load_data(phenotype):
    return (
        load_my_str_results(phenotype),
        load_plink_results(phenotype),
        *load_gwas_catalog(phenotype)
    )

def load_my_str_results(phenotype):
    print(f"Loading my STR results for {phenotype} ... ", end='', flush=True)
    start_time = time.time()
    str_results_fname = f'{ukb}/association/plots/input/{phenotype}/my_str_results.tab'
    my_str_results = df_to_recarray(pd.read_csv(
        str_results_fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        dtype=get_dtypes(str_results_fname)
    ))
    names = list(my_str_results.dtype.names)
    for idx, name in my_results_rename.items():
        names[idx] = name
    for idx, name in my_str_results_rename.items():
        names[idx] = name
    my_str_results.dtype.names = names
    my_str_results['p_val'] = -np.log10(my_str_results['p_val'])
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return my_str_results

def load_my_chr21_height_snp_results():
    print("Loading my chr21 height results ... ", end='', flush=True)
    start_time = time.time()
    snp_results_fname = f'{ukb}/association/plots/input/height/my_imputed_snp_chr21_results.tab'
    my_snp_results = df_to_recarray(pd.read_csv(
        snp_results_fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        dtype=get_dtypes(snp_results_fname)
    ))
    names = list(my_snp_results.dtype.names)
    for idx, name in my_results_rename.items():
        names[idx] = name
    my_snp_results.dtype.names = names
    my_snp_results['p_val'] = -np.log10(my_snp_results['p_val'])
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return my_snp_results

def load_plink_results(phenotype):
    # Load plink SNP results
    print(f"Loading plink SNP results for {phenotype} ... ", end='', flush=True)
    start_time = time.time()
    plink_results = df_to_recarray(pd.read_csv(
        f'{ukb}/association/plots/input/{phenotype}/plink_snp_results.tab',
        sep='\t',
        usecols=(0,1,2,3,4,13,14),
        encoding='UTF-8',
        header=0,
        names=('chr', 'pos', 'id', 'ref', 'alt', 'p_val', 'error')
    ))

    plink_results['p_val'] = -np.log10(plink_results['p_val'])
    plink_results = plink_results[plink_results['p_val'] >= 3]

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--my-plink-comparison", action='store_true')
    parser.add_argument("--phenotype")
    args = parser.parse_args()
    assert bool(args.my_plink_comparison) != bool(args.phenotype)
    if args.phenotype:
        phenotype = args.phentoype
        my_str_results, plink_snp_results, gwas_catalog, gwas_catalog_ids = load_data(
            phenotype
        )
        with open(f'{ukb}/association/results/{phenotype}/my_str/README.txt') as README:
            date_line = next(README)
            my_str_run_date = date_line.split(' ')[2]
    else:
        phenotype = 'height'
        my_snp_results = load_my_chr21_height_snp_results()
        plink_snp_results = load_plink_results(phenotype)
        plink_snp_results = plink_snp_results[plink_snp_results['chr'] == 21]
        gwas_catalog = gwas_catalog_ids = None
        my_str_results = my_str_run_date = None

        with open(f'{ukb}/association/results/{phenotype}/my_str/README.txt') as README:
            date_line = next(README)
            my_snp_run_date = date_line.split(' ')[2]

    with open(f'{ukb}/association/results/{phenotype}/plink_snp/logs/chr21.plink.stdout') as README:
        plink_snp_run_date = ' '.join(README.readlines()[-1].split(' ')[2:])

    with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt') as unit_file:
        unit = next(unit_file).strip()

    make_manhattan_plots(
        phenotype,
        unit,
        my_str_results,
        my_str_run_date,
        plink_snp_results,
        plink_snp_run_date,
        gwas_catalog,
        gwas_catalog_ids,
        my_snp_results,
        my_snp_run_date
    )

if __name__ == "__main__":
    main()

