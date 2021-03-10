#!/usr/bin/env python3

import argparse
import os
import os.path
import time

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
import numpy.ma
import numpy.random
import pandas as pd


import load_and_filter_genotypes

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

def calc_plot_pos(cds: bokeh.models.ColumnDataSource):
    poses = cds.data['pos']
    chroms = cds.data['chr']
    plot_poses = poses.copy()
    for chrom in range(1, 23):
        plot_poses[chroms > chrom] += chr_lens[chrom - 1]
    cds.data['plot_pos'] = plot_poses

def plot_manhattan(plot, source, label, color, size=4):
    return plot.circle(
        'plot_pos',
        'display_p_val',
        legend_label=label,
        source=source,
        color=color,
        size=size,
        muted_alpha=0.1
    )

start_p_val_cap = 30
max_p_val = 350

def make_manhattan_plots(
        phenotype,
        my_str_results,
        my_str_run_date,
        my_snp_results,
        my_snp_run_date,
        plink_snp_results,
        plink_snp_run_name,
        plink_snp_run_date,
        gwas_catalog,
        gwas_catalog_ids):

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

    my_str_data = my_str_results[phenotype]
    my_str_source = bokeh.models.ColumnDataSource(dict())
    my_str_source.data['display_p_val'] = np.minimum(
        my_str_data['p_val'], start_p_val_cap
    )
    for detail_name in my_str_data.dtype.names:
        if detail_name in cols_to_skip:
            continue
        my_str_source.data[detail_name] = my_str_data[detail_name]
    calc_plot_pos(my_str_source)

    my_snp_data = my_snp_results[phenotype]
    my_snp_source = bokeh.models.ColumnDataSource(dict())
    my_snp_source.data['display_p_val'] = np.minimum(
        my_snp_data['p_val'], start_p_val_cap
    )
    for detail_name in ('chr', 'pos', 'alleles', 'p_val'):
        my_snp_source.data[detail_name] = my_snp_data[detail_name]
    calc_plot_pos(my_snp_source)

    plink_snp_data = plink_snp_results[phenotype]
    plink_snp_source = bokeh.models.ColumnDataSource(dict(
        chr=plink_snp_data['chr'],
        pos=plink_snp_data['pos'],
        p_val=plink_snp_data['p_val'],
        display_p_val=np.minimum(plink_snp_data['p_val'], start_p_val_cap),
        alleles=np.char.add(np.char.add(
            plink_snp_data['ref'], ','), plink_snp_data['alt']
        ),
        id=plink_snp_data['id']
    ))
    calc_plot_pos(plink_snp_source)

    catalog_source = bokeh.models.ColumnDataSource(dict(
        chr=gwas_catalog[phenotype][:, 0],
        pos=gwas_catalog[phenotype][:, 1],
        p_val=gwas_catalog[phenotype][:, 2],
        display_p_val=np.minimum(gwas_catalog[phenotype][:, 2], start_p_val_cap),
        rsids=gwas_catalog_ids[phenotype]
    ))
    calc_plot_pos(catalog_source)

    # set up drawing canvas
    plot = bokeh.plotting.figure(
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
    plot.add_tools(wheel_zoom)
    plot.toolbar.active_scroll = wheel_zoom

    pan = bokeh.models.tools.PanTool(dimensions="width")
    plot.add_tools(pan)
    plot.toolbar.active_drag = pan

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
    plot.js_on_event(bokeh.events.MouseEnter, one_tooltip_callback)

    # plot grey backgrounds for even chroms so chroms can be distinguished
    plot.quad(
        top=[max_p_val]*11,
        bottom=[0]*11,
        left=cum_lens[::2],
        right=cum_lens[1::2],
        alpha=0.3,
        color=(220,220,220)
    )

    '''
    line_source = bokeh.models.ColumnDataSource(dict(
        x=[0, chr_lens[1]],
        y=[-np.log10(5e-8)]*2,
    '''
    plot.line(
        x=[0, cum_lens[21]],
        y=[-np.log10(5e-8)]*2,
        line_width=3,
        #line_dash='dashed',
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
    my_snp_manhattan = plot_manhattan(
        plot, my_snp_source, 'SNPs my code', color=(230, 159, 0)
    )
    plink_snp_manhattan = plot_manhattan(
        plot, plink_snp_source, 'SNPs Plink', color=(86, 180, 233)
    )
    my_str_manhattan = plot_manhattan(
        plot, my_str_source, 'STRs my code', color=(204, 121, 167), size=6
    )
    catalog_manhattan = plot_manhattan(
        plot, catalog_source, 'GWAS Catalog Hits', color=(213, 94, 0), size=7
    )

    # add custom tick formatter. See FuncTickFormatter
    # here: https://docs.bokeh.org/en/latest/docs/user_guide/styling.html#userguide-styling-axes-tick-label-formats
    setup_cum_len_js = "var cum_lens = [];"
    for cum_len in cum_lens:
        setup_cum_len_js += f"cum_lens.push({cum_len - chr_lens[0]});"

    plot.xaxis[0].formatter = bokeh.models.formatters.FuncTickFormatter(
        code = setup_cum_len_js + """
        var shift = 0;
        for(var i=0; i<22; i++) {
            if (i > 0) {
                shift = cum_lens[i-1];
            }
            if (tick < cum_lens[i]) {
                return `chr${i}:` + Math.floor((tick - shift)/1000000).toString() + ',' + Math.floor(((tick - shift) % 1000000)/1000).toString() + 'kb';
            }
        }
    """)

    # add hover tooltips for each manhattan plot
    # see https://stackoverflow.com/questions/49282078/multiple-hovertools-for-different-lines-bokeh
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
    plot.add_tools(my_str_hover)

    my_snp_hover = bokeh.models.tools.HoverTool(renderers=[my_snp_manhattan])
    my_snp_hover.tooltips = [
        ('var type', 'SNP'),
        ('alleles:', '@alleles'),
        ('pos', '@pos'),
        ('-log10(p_val) my code', '@p_val')
    ]
    plink_snp_hover = bokeh.models.tools.HoverTool(renderers=[plink_snp_manhattan])
    plink_snp_hover.tooltips = [
        ('var type', 'SNP'),
        ('alleles:', '@alleles'),
        ('pos', '@pos'),
        ('ID', '@id'),
        ('-log10(p_val) Plink', '@p_val')
    ]
    plot.add_tools(plink_snp_hover)

    catalog_hover = bokeh.models.tools.HoverTool(renderers=[catalog_manhattan])
    catalog_hover.tooltips = [
        ('var type', 'GWAS Catalog SNP'),
        ('pos', '@pos'),
        ('rsid:', '@rsids'),
        ('-log10(p_val)', '@p_val')
    ]
    plot.add_tools(catalog_hover)

    plot.toolbar.active_inspect = [my_str_hover, my_snp_hover, plink_snp_hover, catalog_hover]

    height_slider = bokeh.models.Slider(
        start = 8, end=max_p_val, value=start_p_val_cap, step=1, title="p-value cap"
    )
    for source in [my_str_source, my_snp_source, plink_snp_source, catalog_source]:
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
        args=dict(height_slider=height_slider, y_range=plot.y_range),
        code="""
            const new_max = height_slider.value;
            y_range.start = -.25
            y_range.end = new_max + 0.25;
            y_range.change.emit();
        """
    )
    height_slider.js_on_change('value', height_callback)
    layout = bokeh.layouts.column(
        height_slider,
        plot
    )

    plot.add_layout(bokeh.models.Title(
        text=f"My STR code run date: {my_str_run_date}",
        align='right'
    ), 'below')
    plot.add_layout(bokeh.models.Title(
        text=f"My SNP code run date: {my_snp_run_date}",
        align='right'
    ), 'below')
    plot.add_layout(bokeh.models.Title(
        text=f"Plink SNP run name, date: {plink_snp_run_name}, {plink_snp_run_date}",
        align='right'
    ), 'below')

    plot.legend.click_policy="mute"

    html = bokeh.embed.file_html(layout, bokeh.resources.CDN, 'Manhattan plot test')
    with open(f'{ukb}/association/plots/interactive_manhattan.html', 'w') as out_file:
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

def load_data(plink_snp_run_name):
    my_results_rename = {
        0: 'chr',
        3: 'filtered',
        4: 'p_val',
        5: 'coeff_phenotype'
    }

    my_str_results = {}
    print(f"Loading my STR results for {'height'} ... ", end='', flush=True)
    start_time = time.time()
    str_results_fname = f'{ukb}/association/plots/input/{"height"}/my_str_results.tab'
    my_str_results['height'] = df_to_recarray(pd.read_csv(
        str_results_fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        dtype=get_dtypes(str_results_fname)
    ))
    names = list(my_str_results['height'].dtype.names)
    for idx, name in my_results_rename.items():
        names[idx] = name
    my_str_results['height'].dtype.names = names
    my_str_results['height']['p_val'] = -np.log10(my_str_results['height']['p_val'])
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

    my_snp_results = {}
    print(f"Loading my SNP results for {'height'} ... ", end='', flush=True)
    start_time = time.time()
    snp_results_fname=f'{ukb}/association/plots/input/{"height"}/my_imputed_snp_results.tab'
    my_snp_results['height'] = df_to_recarray(pd.read_csv(
        snp_results_fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        dtype=get_dtypes(snp_results_fname)
    ))
    names = list(my_snp_results['height'].dtype.names)
    for idx, name in my_results_rename.items():
        names[idx] = name
    my_snp_results['height'].dtype.names = names
    my_snp_results['height']['p_val'] = -np.log10(my_snp_results['height']['p_val'])
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

    '''
    snp_summary_stats = {}
    snp_ss_description = {}

    print(f"Loading SNP summary stats for height  ... ", end='', flush=True)
    start_time = time.time()
    snp_summary_stats['height'] = np.loadtxt(
        (f"{ukb}/misc_data/snp_summary_stats/height/"
         "Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt"),
        usecols=(0, 1, 8),
        skiprows=1
    )
    snp_ss_description['height'] = 'Ancestry: European, n=700,000'
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

    print(f"Loading SNP summary stats for total bilirubin ... ", end='', flush=True)
    start_time = time.time()
    snp_summary_stats['total_bilirubin'] = np.loadtxt(
        (f"{ukb}/misc_data/snp_summary_stats/bilirubin/"
         "phenocode-TBil_GWAS_in_BBJ_autosome.tsv"),
        usecols=(0, 1, 6),
        skiprows=1,
        delimiter='\t'
    )
    snp_ss_description['total_bilirubin'] = 'Ancestry: Japanese, n=110,000'
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

    for 'height' in snp_summary_stats:
        snp_summary_stats['height'] = prep_data(snp_summary_stats['height'])
    '''

    # Load plink SNP results
    plink_results = {}
    print(f"Loading plink SNP results for {'height'} ... ", end='', flush=True)
    start_time = time.time()
    plink_results['height'] = df_to_recarray(pd.read_csv(
        f'{ukb}/association/runs/{plink_snp_run_name}/results/chr21/plink2.{"height"}_inv_norm_rank.glm.linear.summary',
        sep='\t',
        usecols=(0,1,2,3,4,13,14),
        encoding='UTF-8',
        header=0,
        names=('chr', 'pos', 'id', 'ref', 'alt', 'p_val', 'error')
    ))

    plink_results['height']['p_val'] = -np.log10(plink_results['height']['p_val'])
    plink_results['height'] = plink_results['height'][
        plink_results['height']['p_val'] >= 3
    ]

    assert np.all(plink_results['height']['error'] == '.')

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

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

    for phenotype in known_assocs:
        zero_idx = known_assocs[phenotype][:, 2] == 0
        known_assocs[phenotype][~zero_idx, 2] = -np.log10(known_assocs[phenotype][~zero_idx, 2])
        known_assocs[phenotype][zero_idx, 2] = 50.12345 #choose an arbitrary number

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

    return my_str_results, my_snp_results, plink_results, known_assocs, known_assoc_ids


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("plink_snp_run_name")
    args = parser.parse_args()

    my_str_results, my_snp_results, plink_snp_results, gwas_catalog, gwas_catalog_ids = load_data(
        args.plink_snp_run_name
    )
    print(my_str_results['height'].shape,
          my_snp_results['height'].shape, plink_snp_results['height'].shape, gwas_catalog['height'].shape)

    with open(f'{ukb}/association/results/height/my_str/README.txt') as README:
        next(README)
        date_line = next(README)
        my_str_run_date = date_line.split(' ')[2]

    with open(f'{ukb}/association/results/height/my_imputed_snp/README.txt') as README:
        next(README)
        date_line = next(README)
        my_snp_run_date = date_line.split(' ')[2]

    with open(f'{ukb}/association/runs/{args.plink_snp_run_name}/README') as README:
        next(README)
        date_line = next(README)
        plink_snp_run_date = date_line.split(' ')[2]

    make_manhattan_plots(
        'height',
        my_str_results,
        my_str_run_date,
        my_snp_results,
        my_snp_run_date,
        plink_snp_results,
        args.plink_snp_run_name,
        plink_snp_run_date,
        gwas_catalog,
        gwas_catalog_ids
    )

if __name__ == "__main__":
    main()

