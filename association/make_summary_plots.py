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

def plot_manhattan(ax, data, label, color):
    # Trim to chroms 2, 3
    data = data[np.logical_or(data[:, 0] == 2, data[:, 0] == 3), :]
    xs = data[:, 1]
    xs[data[:, 0] == 3] += 243199373 #offset by len of chr 2
    ax.scatter(
        xs,
        data[:, 2],
        color=colors[color],
        label=label,
        marker='.'
    )

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


def make_manhattan_plots(phenotypes, results, snp_summary_stats, run_name):
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
        ax.set_xlabel(f'Position (bp count starts on chr2, continues into chr3)')

        max_len = 243199373 + 198022430 # len chr2 + chr3
        ax.plot([0, max_len], [GENOME_WIDE_TRANS_SIG] * 2, label='p >= 5e-8',
                 color=colors['red'])

        if phenotype in snp_summary_stats:
            plot_manhattan(ax, snp_summary_stats[phenotype], 'SNPs', 'royalblue')

        plot_manhattan(ax, results[phenotype], 'STRs', 'mediumvioletred')

        ax.legend()

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

    make_qq_plots(phenotypes, results, snp_summary_stats, args.run_name)
    make_manhattan_plots(phenotypes, results, snp_summary_stats,
                         args.run_name)

if __name__ == "__main__":
    main()

