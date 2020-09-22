import argparse
import os
import time

from matplotlib.colors import CSS4_COLORS as colors
import matplotlib.colors
import matplotlib.image
import matplotlib.patches
import matplotlib.path
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import numpy as np
import numpy.ma
import numpy.random

import scipy.ndimage

import load_and_filter_genotypes

ukb = os.environ['UKB']
GENOME_WIDE_TRANS_SIG = -np.log10(5e-8)
HEIGHT_CUTOFF = 20

# assume data.shape is (loci, 3)
# where columns are chrom, pos, pval
# will return sorted from least to most significant
# with values in the transformed axis (-np.log10(pval))
def prep_data(data):
    data = data[~np.isnan(data[:, 2]), :]
    pvals = data[:, 2]
    pvals[pvals == 0] = min(pvals[pvals != 0])
    pvals[:] = -np.log10(pvals)
    return data[np.argsort(data[:, 2])]


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


def plot_manhattan(ax, data, label, color_pair, chr_lens, height_cutoff):
    if height_cutoff:
        data = data.copy()
        data[data[:, 2] > HEIGHT_CUTOFF, 2] = HEIGHT_CUTOFF
    first = True
    for chrom in range(1, 23):
        if not first:
            label = None
        chrom_data = data[data[:, 0] == chrom, :]
        xs = chrom_data[:, 1]
        xs += np.sum(chr_lens[0:(chrom-1)])
        ax.scatter(
            xs,
            chrom_data[:, 2],
            color=colors[color_pair[chrom % 2]],
            label=label,
            marker='.'
        )
        ax.margins(0, 0.05)
        first = False


def make_genome_manhattan_plots(phenotypes, results, snp_summary_stats,
                                snp_ss_description, known_assocs, run_name,
                                height_cutoff):
    # plot Manhattan plots
    nphenotypes = len(phenotypes)
    fig, axs = plt.subplots(nphenotypes, figsize=(300, 5*nphenotypes))
    fig.suptitle('Manhattan plots')
    fig.tight_layout(pad=5.0)
    chr_lens = np.genfromtxt(
        f'{ukb}/misc_data/genome/chr_lens.txt',
        skip_header=1,
        usecols=(1),
        dtype=int
    )
    max_len = np.sum(chr_lens)

    plt_count = 0
    for phenotype in phenotypes:
        if len(phenotypes) > 1:
            ax = axs[plt_count]
        else:
            ax = axs
        ax.set_title(f'{phenotype}')
        ax.set_ylabel('-log_10(pval)')
        ax.set_xlabel('Position (1e7 bp)')

        ax.plot([0, max_len], [GENOME_WIDE_TRANS_SIG] * 2, label='p = 5e-8',
                 color=colors['red'])

        if phenotype in snp_summary_stats:
            desc = f'SNPs ({snp_ss_description[phenotype]})'
            plot_manhattan(ax, snp_summary_stats[phenotype], desc,
                           ('royalblue', 'cornflowerblue'), chr_lens,
                           height_cutoff)

        plot_manhattan(ax, results[phenotype],
                       'STRs (Ancestry: European, n=500,000)',
                       ('mediumvioletred', 'deeppink'), chr_lens,
                       height_cutoff)

        if phenotype in known_assocs:
            plot_manhattan(ax, known_assocs[phenotype], 'GWAS Catalog Hits',
                           ('goldenrod', 'gold'), chr_lens, height_cutoff)


        ax.legend()

        ticks = []
        tick_labels = []
        for chrom, length in enumerate(chr_lens):
            chrom += 1
            new_ticks = np.arange(0, length, int(5e7))
            new_ticks += np.sum(chr_lens[0:(chrom-1)])
            ticks.extend(new_ticks)
            tick_labels.append(f'chr{chrom}:0')
            tick_labels.extend(str(int(tick)) for tick in np.arange(5, length/1e7, 5))
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels)

        plt_count += 1

    cutoff_txt = ''
    if not height_cutoff:
        cutoff_txt = '_no_cutoff'
    plt.savefig(f'{ukb}/association/runs/{run_name}/manhattan_plot_all{cutoff_txt}.png')


def make_region_plots(results, snp_summary_stats,
                     snp_ss_description, known_assocs, run_name,
                     imputation_run_name, filtering_run_name, region):
    chrom, locus, phenotype = region.split(":")
    chrom = int(chrom)
    if '-' in locus:
        start, end = (int(x) for x in locus.split('-'))
    else:
        pos = int(locus)
        start, end = pos - int(5e5), pos + int(5e5)

    os.makedirs(f'{ukb}/association/runs/{run_name}/regional_plots', exist_ok=True)

    fname = f'{ukb}/sample_qc/runs/{filtering_run_name}/combined_unrelated.sample'
    samplelist = np.genfromtxt(
        fname,
        skip_header=1,
        usecols=[0],
        dtype='U7'
    ).reshape(-1)
    samplelist = np.char.add(np.char.add(samplelist, '_'), samplelist)

    gts_per_str, str_poses, str_samples  = load_and_filter_genotypes.load_all_haplotypes(
        load_and_filter_genotypes.filtered_strs(
            imputation_run_name,
            f'{chrom}:{start}-{end}'
        ),
        samplelist
    )
    gts_per_snp, snp_poses, snp_samples = load_and_filter_genotypes.load_all_haplotypes(
        load_and_filter_genotypes.filtered_snps(f'{chrom}:{start}-{end}'),
        samplelist
    )
    assert set(str_samples) == set(snp_samples)

    for plot_num in range(3):
        if plot_num == 0:
            print("Making STR heatmap plot ... ", flush=True)
        elif plot_num == 1:
            print("Making SNP heatmap plot ... ", flush=True)
        elif plot_num == 2:
            print("Making SNP-STR heatmap plot ... ", flush=True)
        else:
            raise ValueError()

        fig_height = 14
        fig_width = 14
        width_sides = 0.1
        height_sides = width_sides * fig_width/fig_height
        absolute_height_spacing = 0.15 * 8
        height_spacing = absolute_height_spacing/fig_height
        absolute_width_spacing = 0.8
        width_spacing = absolute_width_spacing/fig_width
        cbar_width_abs = 0.25
        cbar_width = cbar_width_abs/fig_width
        fig = plt.figure(figsize = (fig_width, fig_height))
        manhattan_ax = fig.add_axes([
            width_sides,
            fig_width/fig_height*(1 - 2*width_sides - width_spacing - cbar_width)/np.sqrt(2) + height_spacing,
            1 - 2*width_sides - width_spacing - cbar_width,
            1 - (fig_width/fig_height*(1 - 2*width_sides - width_spacing - cbar_width)/np.sqrt(2) + height_spacing) - fig_width/fig_height*width_sides,
        ])
        heatmap_ax = fig.add_axes([
            width_sides,
            height_sides,
            1 - 2*width_sides - width_spacing - cbar_width,
            fig_width/fig_height*(1 - 2*width_sides - width_spacing - cbar_width)/np.sqrt(2)
        ])
        cbar_ax = fig.add_axes([
            1 - width_sides - cbar_width,
            height_sides,
            cbar_width,
            fig_width/fig_height*(1 - 2*width_sides - width_spacing - cbar_width)/np.sqrt(2)
        ])
        
        fig.suptitle(f'Associations with {phenotype}', fontsize='large')

        print("Drawing Manhattan plot ... ", end='', flush=True)
        start_time = time.time()
        make_region_manhattan_plot(
            manhattan_ax,
            results,
            snp_summary_stats,
            snp_ss_description,
            known_assocs,
            chrom,
            start,
            end,
            phenotype
        )
        print(f"done ({time.time() - start_time:.2e}s)", flush=True)

        if plot_num == 0:
            gts_per_locus = gts_per_str
            poses = str_poses
            fname = 'str'
            cbar_name = 'STRs'
        elif plot_num == 1:
            gts_per_locus = gts_per_snp
            poses = snp_poses
            fname = 'snp'
            cbar_name = 'SNPs'
        elif plot_num == 2:
            def double_array(arr):
                '''1D array arr'''
                arr_2d = arr.reshape(-1, 1)
                return np.concatenate((arr_2d, arr_2d), axis=1).reshape(-1)

            str_sample_sort = np.argsort(str_samples)
            snp_sample_sort = np.argsort(snp_samples)
            gts_per_locus = numpy.ma.concatenate(
                (
                    gts_per_str[double_array(str_sample_sort), :],
                    gts_per_snp[double_array(snp_sample_sort), :]
                ),
                axis = 1
            )

            poses = np.concatenate((str_poses, snp_poses))
            fname = 'combined'
            cbar_name = 'loci (STRs and SNPs)'
        else:
            raise ValueError()
        make_ld_heatmap(
            heatmap_ax,
            cbar_ax,
            gts_per_locus,
            poses,
            start,
            end,
            cbar_name
        )

        plt.savefig(f'{ukb}/association/runs/{run_name}/regional_plots/{phenotype}_{chrom}_{start}_{end}_{fname}_heatmap.png')

def make_region_manhattan_plot(ax, results, snp_summary_stats,
                         snp_ss_description, known_assocs,
                         chrom, start, end, phenotype):
    chr_lens = np.genfromtxt(
        f'{ukb}/misc_data/genome/chr_lens.txt',
        skip_header=1,
        usecols=(1),
        dtype=int
    )

    def subset(arr, chrom, start, end):
        arr = arr[arr[:, 0] == chrom, :]
        arr = arr[arr[:, 1] >= start, :]
        return arr[arr[:, 1] <= end, :]

    ax.set_ylabel('-log_10(pval)')
    ax.set_xlabel('Position')

    offset = np.sum(chr_lens[0:(chrom-1)])
    ax.plot([start+offset, end+offset], [GENOME_WIDE_TRANS_SIG] * 2, label='p = 5e-8',
             color=colors['red'])

    if phenotype in snp_summary_stats:
        desc = f'SNPs ({snp_ss_description[phenotype]})'
        plot_manhattan(ax,
                       subset(snp_summary_stats[phenotype], chrom, start, end),
                       desc,
                       ('royalblue', 'cornflowerblue'),
                       chr_lens, False)

    plot_manhattan(ax,
                   subset(results[phenotype], chrom, start, end),
                   'STRs (Ancestry: European, n=500,000)',
                   ('mediumvioletred', 'deeppink'),
                   chr_lens, False)

    if phenotype in known_assocs:
        plot_manhattan(ax,
                       subset(known_assocs[phenotype], chrom, start, end),
                       'GWAS Catalog Hits',
                       ('goldenrod', 'gold'),
                       chr_lens, False)
    ax.legend()


def make_ld_heatmap(ax, cbar_ax, gts_per_locus, poses, start, end, cbar_name):
    rand = numpy.random.default_rng(seed=13)
    gts_per_locus = gts_per_locus[
        rand.choice(gts_per_locus.shape[0], size=int(5e4), replace=False),
        :
    ].copy() # random subset to make computations not take forever
    # perhaps copying this into a new array which is not a view will
    # make the correlation below faster

    print("Correlating and sorting loci ... ", end='', flush=True)
    start_time = time.time()
    corrcoefs = np.square(numpy.ma.corrcoef(gts_per_locus, rowvar=False))

    sort_order = np.argsort(poses)
    poses = poses[sort_order]
    corrcoefs = corrcoefs[sort_order, :]
    corrcoefs = corrcoefs[:, sort_order]
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

    print("Drawing LD heatmap ... ", end='', flush=True)
    start_time = time.time()

    mesh1D = np.linspace(start, end, 10000, endpoint=True)
    meshx, meshy = np.meshgrid(mesh1D, mesh1D)
    def find_closest(mesh):
        nearest_above = np.searchsorted(poses, mesh)
        nearest_below = nearest_above - 1
        nearest_above[nearest_above == len(poses)] = len(poses) - 1
        nearest_below[nearest_below < 0] = 0
        pos_above = poses[nearest_above]
        pos_below = poses[nearest_below]
        nearest_idx = nearest_above
        below_closer = np.abs(pos_below - mesh) < np.abs(pos_above - mesh)
        nearest_idx[below_closer] = nearest_below[below_closer]
        return nearest_idx
    closestx = find_closest(meshx)
    closesty = find_closest(meshy)
    z = corrcoefs[closestx, closesty]
    z = z.reshape(meshx.shape)
    z = scipy.ndimage.rotate(z, 45)

    im = ax.imshow(
        z,
        extent=(start,end,start,end),
        vmin=0,
        vmax=1,
        cmap='YlOrRd',
        norm=matplotlib.colors.TwoSlopeNorm(vmin=0, vcenter=0.05, vmax=1)
    )
    midpoint = (start+end)/2
    im.set_clip_path(matplotlib.patches.PathPatch(
        matplotlib.path.Path(
            [[start, midpoint], [midpoint, start], [end, midpoint]]
        ), transform=ax.transData
    ))

    ax.set_xlim(start, end)
    ax.set_ylim(start, (start+end)/2)
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')

    cbar = ax.figure.colorbar(
        im,
        cax=cbar_ax
    )
    cbar.ax.set_ylabel(f"Correlation between {cbar_name}", rotation=-90, va="bottom")\

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_name")
    parser.add_argument("imputation_run_name")
    parser.add_argument("filtering_run_name")
    parser.add_argument(
        "--regions",
        help=("comma separated list of regions to make summary plots for. "
              "Each region is of the form chrom:start-end:phenotype or "
              "chrom:pos:phenotype in which case a 1Mb region around the "
              "position is shown ")
    )
    args = parser.parse_args()
    if args.regions:
        regions = args.regions.split(',')
    else:
        regions = None

    phenotypes = ['height', 'total_bilirubin']

    results = {}
    for phenotype in phenotypes:
        print(f"Loading results for {phenotype} ... ", end='', flush=True)
        start_time = time.time()
        results[phenotype] = np.loadtxt(
            f'{ukb}/association/runs/{args.run_name}/results/{phenotype}.txt',
            usecols=[0, 1, 4],
            skiprows=1
        )
        print(f"done ({time.time() - start_time:.2e}s)", flush=True)

    #TODO fix nan alleles!
    for phenotype in phenotypes:
        results[phenotype] = prep_data(results[phenotype])

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
    
    # Load the NHGRI-EBI GWAS catalog
    print(f"Loading GWAS catalog results ... ", end='', flush=True)
    start_time = time.time()
    catalog = np.loadtxt(
        f'{ukb}/misc_data/snp_summary_stats/catalog/catalog_hg19.tsv',
        usecols=[11, 12, 27],
        delimiter='\t',
        skiprows=1,
        dtype=object
    )
    # omit results that are not mapped to a recognizable chromosome
    catalog_names = np.loadtxt(
        f'{ukb}/misc_data/snp_summary_stats/catalog/catalog_hg19.tsv',
        usecols=[7],
        delimiter='\t',
        skiprows=1,
        dtype=object
    )
    filter_weird_chroms = np.isin(catalog[:, 0],
                                  list(str(chrom) for chrom in range(1, 23)))
    catalog_names = catalog_names[filter_weird_chroms]
    catalog = catalog[filter_weird_chroms, :]
    catalog = catalog.astype(float)
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

    known_assocs = {}
    known_assocs['height'] = catalog[np.equal(catalog_names, 'Height'), :]
    bil_catalog_trait_names = [
        'Total bilirubin levels',
        'Bilirubin levels'
    ]
    known_assocs['total_bilirubin'] = \
            catalog[np.isin(catalog_names, bil_catalog_trait_names), :]
    for phenotype in known_assocs:
        known_assocs[phenotype] = prep_data(known_assocs[phenotype])

    if not regions:
        for cutoff in {True, False}:
            make_qq_plots(phenotypes, results, snp_summary_stats,
                          args.run_name, cutoff)
            make_genome_manhattan_plots(
                phenotypes,
                results,
                snp_summary_stats,
                snp_ss_description,
                known_assocs,
                args.run_name,
                cutoff
            )
    else:
        for region in regions:
            make_region_plots(
                results,
                snp_summary_stats,
                snp_ss_description,
                known_assocs,
                args.run_name,
                args.imputation_run_name,
                args.filtering_run_name,
                region
            )

if __name__ == "__main__":
    main()

