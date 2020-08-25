import argparse
import os

from matplotlib.colors import CSS4_COLORS as colors
import matplotlib.pyplot as plt
import numpy as np


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
    pvals[pvals > HEIGHT_CUTOFF] = HEIGHT_CUTOFF
    return data[np.argsort(data[:, 2])]

# gets an array of x-axis vals of the corresponding length
def qq_x_axis(length):
    uniform = (np.arange(length) + 0.5)/length
    trans_uniform = -np.log10(uniform)[::-1]
    return trans_uniform

def plot_qq(ax, prepped_pvals, label, color):
    xs = qq_x_axis(len(prepped_pvals))
    ax.scatter(
        xs,
        prepped_pvals,
        color=colors[color],
        marker='.',
        label=label
    )

def plot_manhattan(ax, data, label, color_pair):
    # Trim to chroms 2, 3
    first = True
    for chrom in 2, 3:
        if not first:
            label = None
        chrom_data = data[data[:, 0] == chrom, :]
        xs = chrom_data[:, 1]
        if chrom == 3:
            xs += 243199373 #offset by len of chr 2
        ax.scatter(
            xs,
            chrom_data[:, 2],
            color=colors[color_pair[chrom % 2]],
            label=label,
            marker='.'#,
            #alpha=0.5
        )
        first = False

def make_qq_plots(phenotypes, results, snp_summary_stats, run_name):
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
        ax.set_xlabel(f'-log_10(Uniform dist) (total STR loci = {nresults})')

        if phenotype in snp_summary_stats:
            max_len = max(snp_summary_stats[phenotype].shape[0], nresults)
            plot_qq(ax, snp_summary_stats[phenotype][:, 2],
                    'SNPs genome wide', 'royalblue')
        else:
            max_len = nresults

        plot_qq(ax, results[phenotype][:, 2], 'STRs on chr2,3',
                'mediumvioletred')
        plot_qq(ax, qq_x_axis(max_len),
                'y=x, expectation under the null', 'gray')

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

    plt.savefig(f'{ukb}/association/runs/{run_name}/qq_plot.png')


def make_manhattan_plots(phenotypes, results, snp_summary_stats, known_assocs, run_name):
    # plot Manhattan plots
    nphenotypes = len(phenotypes)
    fig, axs = plt.subplots(nphenotypes, figsize=(30, 5*nphenotypes))
    fig.suptitle('Manhattan plots (chr2,3)')
    fig.tight_layout(pad=5.0)
    plt_count = 0
    for phenotype in phenotypes:
        if len(phenotypes) > 1:
            ax = axs[plt_count]
        else:
            ax = axs
        ax.set_title(f'{phenotype}')
        ax.set_ylabel('-log_10(pval)')
        ax.set_xlabel(f'Position (1e7 bp)')

        max_len = 243199373 + 198022430 # len chr2 + chr3
        ax.plot([0, max_len], [GENOME_WIDE_TRANS_SIG] * 2, label='p = 5e-8',
                 color=colors['red'])

        if phenotype in snp_summary_stats:
            plot_manhattan(ax, snp_summary_stats[phenotype], 'SNPs',
                           ('royalblue', 'cornflowerblue'))

        plot_manhattan(ax, results[phenotype], 'STRs',
                       ('mediumvioletred', 'deeppink'))

        if phenotype in known_assocs:
            plot_manhattan(ax, known_assocs[phenotype], 'GWAS Catalog Hits',
                           ('goldenrod', 'gold'))


        ax.legend()

        ticks = []
        tick_labels = []
        for chrom, length in (2, 243199373), (3, 198022430):
            new_ticks = np.arange(0, length, int(5e7))
            if chrom == 3:
                new_ticks += 243199373
            ticks.extend(new_ticks)
            tick_labels.append(f'chr{chrom}:0')
            tick_labels.extend(str(int(tick)) for tick in np.arange(5, length/1e7, 5))
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels)

        plt_count += 1

    plt.savefig(f'{ukb}/association/runs/{run_name}/manhattan_plot.png')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_name")
    args = parser.parse_args()

    phenotypes = ['height', 'total_bilirubin']

    results = {}
    for phenotype in phenotypes:
        results[phenotype] = np.loadtxt(
            f'{ukb}/association/runs/{args.run_name}/results/{phenotype}.txt',
            usecols=[0, 1, 4],
            skiprows=1
        )
    #TODO fix nan alleles!
    for phenotype in phenotypes:
        results[phenotype] = prep_data(results[phenotype])

    snp_summary_stats = {}
    snp_summary_stats['height'] = np.loadtxt(
        (f"{ukb}/misc_data/snp_summary_stats/height/"
         "Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt"),
        usecols=(0, 1, 8),
        skiprows=1
    )
    snp_summary_stats['total_bilirubin'] = np.loadtxt(
        (f"{ukb}/misc_data/snp_summary_stats/bilirubin/"
         "phenocode-TBil_GWAS_in_BBJ_autosome.tsv"),
        usecols=(0, 1, 6),
        skiprows=1,
        delimiter='\t'
    )
    for phenotype in snp_summary_stats:
        snp_summary_stats[phenotype] = prep_data(snp_summary_stats[phenotype])

    # Load the NHGRI-EBI GWAS catalog
    catalog = np.loadtxt(
        f'{ukb}/misc_data/snp_summary_stats/catalog/catalog.tsv',
        usecols=[11, 12, 27],
        delimiter='\t',
        skiprows=1,
        dtype=object
    )
    # omit results that are not mapped to a recognizable chromosome
    catalog_names = np.loadtxt(
        f'{ukb}/misc_data/snp_summary_stats/catalog/catalog.tsv',
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


    # make_qq_plots(phenotypes, results, snp_summary_stats, args.run_name)
    make_manhattan_plots(phenotypes, results, snp_summary_stats, known_assocs,
                        args.run_name)

if __name__ == "__main__":
    main()

