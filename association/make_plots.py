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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_name")
    args = parser.parse_args()

    phenotypes = ['height', 'total_bilirubin']
    #phenotypes = ['height']

    results = {}
    results['height'] = np.loadtxt(
        f'{ukb}/association/runs/{args.run_name}/results.txt',
        usecols=[0, 1, 4],
        skiprows=1
    )
    results['total_bilirubin'] = np.loadtxt(
        f'{ukb}/association/runs/{args.run_name}/results.txt',
        usecols=[0, 1, 6],
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
    for phenotype in snp_summary_stats:
        snp_summary_stats[phenotype] = prep_data(snp_summary_stats[phenotype])

    # plot QQ plots
    n_phenotypes = len(phenotypes)
    fig, axs = plt.subplots(n_phenotypes, figsize=(10, 5*n_phenotypes))
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

    plt.savefig(f'{ukb}/association/runs/{args.run_name}/qq_plot.png')


if __name__ == "__main__":
    main()

