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


import load_and_filter_genotypes

ukb = os.environ['UKB']
GENOME_WIDE_TRANS_SIG = -np.log10(5e-8)
HEIGHT_CUTOFF = 20

# gets an array of x-axis vals of the corresponding length
def qq_x_axis(length):
    uniform = (np.arange(length) + 0.5)/length
    trans_uniform = -np.log10(uniform)[::-1]
    return trans_uniform


def plot_qq(ax, prepped_pvals, label, color, height_cutoff):
    if height_cutoff:
        prepped_pvals = prepped_pvals.copy()
        prepped_pvals[prepped_pvals > HEIGHT_CUTOFF] = HEIGHT_CUTOFF
    xs = qq_x_axis(len(prepped_pvals))
    ax.scatter(
        xs,
        prepped_pvals,
        color=colors[color],
        marker='.',
        label=label
    )


def make_qq_plots(phenotypes, results, snp_summary_stats, run_name,
                  height_cutoff):
    # plot QQ plots
    nphenotypes = len(phenotypes)
    fig, axs = plt.subplots(nphenotypes, figsize=(10, 5*nphenotypes))
    fig.suptitle('QQ plots')
    fig.tight_layout(pad=5.0)
    plt_count = 0
    for phenotype in phenotypes:
        nresults = results[phenotype].shape[0]
        if len(phenotypes) > 1:
            ax = axs[plt_count]
        else:
            ax = axs
        ax.set_title(f'{phenotype}')
        ax.set_ylabel('-log_10(pval)')
        ax.set_xlabel(f'-log_10(Uniform dist) (# STR loci = {nresults})')

        if phenotype in snp_summary_stats:
            max_len = max(snp_summary_stats[phenotype].shape[0], nresults)
            plot_qq(ax, snp_summary_stats[phenotype][:, 2],
                    'SNPs', 'royalblue',
                    height_cutoff)
        else:
            max_len = nresults

        plot_qq(ax, results[phenotype][:, 2], 'STRs',
                'mediumvioletred', height_cutoff)
        plot_qq(ax, qq_x_axis(max_len),
                'y=x, expectation under the null', 'gray',
                height_cutoff)

        crosses_gws = np.argmax(results[phenotype][:, 2]
                                > GENOME_WIDE_TRANS_SIG)
        xy = (qq_x_axis(nresults)[crosses_gws],
              results[phenotype][crosses_gws, 2])
        n_sig_loci = np.sum(results[phenotype][:, 2]
                            > GENOME_WIDE_TRANS_SIG)
        ax.annotate(
            f"{n_sig_loci} loci with p <= 5e-8",
            xy=xy,
            xytext=(xy[0]-1.5, xy[1]+5),
            arrowprops=dict(facecolor='black', shrink=0.05),
        )

        ax.legend()

        plt_count += 1

    cutoff_txt = ''
    if not height_cutoff:
        cutoff_txt = '_no_cutoff'
    plt.savefig(f'{ukb}/association/runs/{run_name}/qq_plot{cutoff_txt}.png')


def plot_manhattan(plot, source, label, chr_lens):
    assert np.all(source.data['chr'] == 21)
    plot.circle(
        'pos',
        'display_p_val',
        legend_label=label,
        source=source
    )

def make_manhattan_plots(phenotypes,
                         my_snp_results,
                         my_snp_run_name,
                         my_snp_run_date,
                         gwas_catalog,
                         gwas_catalog_ids):
    # plot Manhattan plots
    chr_lens = np.genfromtxt(
        f'{ukb}/misc_data/genome/chr_lens.txt',
        skip_header=1,
        usecols=(1),
        dtype=int
    )

    for phenotype in phenotypes:
        print(f"Plotting phenotype {phenotype} ... ")
        my_snp_data = my_snp_results[phenotype]
        # reformat data for plotting
        start_p_val_cap = 30
        my_snp_source = bokeh.models.ColumnDataSource(dict(
            chr=my_snp_data['chr'],
            pos=my_snp_data['pos'],
            p_val=my_snp_data['p_val'],
            display_p_val=np.minimum(my_snp_data['p_val'], start_p_val_cap),
            alleles=my_snp_data['alleles']
        ))
        snp_detail_names = my_snp_data.dtype.names[7:]
        for detail_name in snp_detail_names:
            my_snp_source.data[detail_name] = my_snp_data[detail_name]

        catalog_source = bokeh.models.ColumnDataSource(dict(
            chr=gwas_catalog[phenotype][:, 0],
            pos=gwas_catalog[phenotype][:, 1],
            p_val=gwas_catalog[phenotype][:, 2],
            rsids=gwas_catalog_ids[phenotype]
        ))

        min_x = min(my_snp_data['pos'])
        max_x = max(my_snp_data['pos'])

        x_diff = max_x - min_x
        margin = 0.15

        # set up drawing canvas
        plot = bokeh.plotting.figure(
            width=1200,
            height=500,
            title=f'{phenotype} Manhattan plots'.capitalize(),
            x_axis_label='chr21 position',
            y_axis_label='-log10(p-value)',
            x_range = bokeh.models.DataRange1d(
                min_x,
                max_x,
                bounds = (min_x - margin*x_diff, max_x + margin*x_diff)
            ),
            y_range = bokeh.models.DataRange1d(
                0,
                start_p_val_cap,
                bounds = None
            ),
            tools='xzoom_in,xzoom_out,save'
        )

        # add custom tools
        wheel_zoom = bokeh.models.tools.WheelZoomTool(dimensions="width")
        plot.add_tools(wheel_zoom)
        plot.toolbar.active_scroll = wheel_zoom

        pan = bokeh.models.tools.PanTool(dimensions="width")
        plot.add_tools(pan)
        plot.toolbar.active_drag = pan

        hover = bokeh.models.tools.HoverTool()
        hover.tooltips = [
            ('var type', 'SNP'),
            ('alleles:', '@alleles'),
            ('pos', '@pos'),
            ('p_val my code', '@p_val')
        ]
        for detail_name in snp_detail_names:
            hover.tooltips.append((detail_name, f'@{detail_name}'))

        plot.add_tools(hover)
        plot.toolbar.active_inspect = [hover]

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

        plot.line(
            x=[0, chr_lens[21]],
            y=[-np.log(5e-8)]*2,
            line_width=3,
            line_dash='dashed',
            color='red',
            legend_label='GWAS p-value threshold'
        )
        plot_manhattan(
            plot, my_snp_source, 'SNPs my code', chr_lens
        )
        plot_manhattan(


        height_slider = bokeh.models.Slider(
            start = 8, end=350, value=start_p_val_cap, step=1, title="p-value cap"
        )
        for source in [my_snp_source]:
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
                y_range.end = new_max + 1;
                y_range.change.emit();
            """
        )
        height_slider.js_on_change('value', height_callback)
        layout = bokeh.layouts.column(
            height_slider,
            plot
        )

        plot.add_layout(bokeh.models.Title(
            text=f"My code run name: {my_snp_run_name}",
            align='right'
        ), 'below')
        plot.add_layout(bokeh.models.Title(
            text=f"My code run date: {my_snp_run_date}",
            align='right'
        ), 'below')

        html = bokeh.embed.file_html(layout, bokeh.resources.CDN, 'Manhattan plot test')
        with open(f'{ukb}/association/plots/me_manhattan_{phenotype}.html', 'w') as out_file:
            out_file.write(html)

    return

def load_data(phenotypes, my_snp_run_name):
    my_snp_results = {}
    for phenotype in phenotypes:
        print(f"Loading results for {phenotype} ... ", end='', flush=True)
        start_time = time.time()
        my_snp_results[phenotype] = np.genfromtxt(
            f'{ukb}/association/runs/{my_snp_run_name}/results/{phenotype}.tab',
            names=True,
            delimiter='\t',
            dtype=None,
            encoding='UTF-8'
        )
        names = list(my_snp_results[phenotype].dtype.names)
        new_names = ('chr', 'pos', 'alleles', 'filtered', 'p_val',
                     'coeff_height', 'coeff_intercept')
        for idx, new_name in enumerate(new_names):
            names[idx] = new_name
        my_snp_results[phenotype].dtype.names = names
        my_snp_results[phenotype]['p_val'] = -np.log10(my_snp_results[phenotype]['p_val'])
        my_snp_results[phenotype] = my_snp_results[phenotype][
            my_snp_results[phenotype]['p_val'] >= 3
        ]
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

    for phenotype in snp_summary_stats:
        snp_summary_stats[phenotype] = prep_data(snp_summary_stats[phenotype])
    '''
    
    # Load the NHGRI-EBI GWAS catalog
    print(f"Loading GWAS catalog results ... ", end='', flush=True)
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
    ).astype('U') # cols 7, 34 describe the phenotype, 21 is the snp rsid
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
            catalog[catalog_names[:, 1] == 'bilirubin measurement', :]
    known_assoc_ids['total_bilirubin'] = \
            catalog_names[catalog_names[:, 1] == 'bilirubin measurement', 1]

    for phenotype in known_assocs:
        known_assocs[phenotype][:, 2] = -np.log10(known_assocs[phenotype][:, 2])
        known_assocs[phenotype] = known_assocs[phenotype][
            known_assocs[phenotype][:, 0] == '21',
            :
        ]
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

    return my_snp_results, known_assocs, known_assoc_ids


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("my_snp_run_name")
    args = parser.parse_args()
    phenotypes = ['height', 'total_bilirubin']

    my_snp_results, gwas_catalog, gwas_catalog_ids = load_data(phenotypes, args.my_snp_run_name)

    with open(f'{ukb}/association/runs/{args.my_snp_run_name}/README') as README:
        next(README)
        date_line = next(README)
        my_snp_run_date = date_line.split(' ')[2]

    make_manhattan_plots(
        phenotypes,
        my_snp_results,
        args.my_snp_run_name,
        my_snp_run_date,
        gwas_catalog,
        gwas_catalog_ids
    )

if __name__ == "__main__":
    main()

