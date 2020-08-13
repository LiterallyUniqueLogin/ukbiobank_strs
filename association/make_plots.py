import argparse
import os

from matplotlib.colors import CSS4_COLORS as colors
import matplotlib.pyplot as plt
import numpy as np


ukb = os.environ['UKB']
GENOME_WIDE_SIG = 5e-8
HEIGHT_CUTOFF = 20

def plot_height_snp_qq(ax):
    snp_summary_stats = np.loadtxt(
        (f"{ukb}/misc_data/snp_summary_stats/height/"
         "Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt"),
        skiprows=1,
        usecols=(8)
    ).reshape(-1)
    snp_summary_stats[snp_summary_stats == 0] = \
            min(snp_summary_stats[snp_summary_stats != 0])
    snp_summary_stats.sort()
    trans_data = -np.log10(snp_summary_stats)
    trans_data[trans_data > HEIGHT_CUTOFF] = HEIGHT_CUTOFF
    uniform = (np.arange(len(trans_data)) + 0.5)/len(trans_data)
    trans_uniform = -np.log10(uniform)
    ax.scatter(
        trans_uniform,
        trans_data,
        color=colors['royalblue'],
        marker='.',
        label="SNPs genome wide"
    )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_name")
    args = parser.parse_args()

    results = np.loadtxt(
        f'{ukb}/association/runs/{args.run_name}/results.txt',
        skiprows=1,
        usecols=(0, 1, 4, 6, 8, 10)
    )
    pvals = {}
    pvals['height'] = results[:, 2]
    pvals['total_bilirubin'] = results[:, 3]
    """
    pvals['direct_bilirubin'] = results[:, 4]
    pvals['indirect_bilirubin'] = results[:, 5]
    """


    nresults = results.shape[0]
    uniform = (np.arange(nresults) + 0.5)/nresults
    trans_uniform = -np.log10(uniform)

    fig, axs = plt.subplots(2, figsize=(10, 10))
    fig.suptitle('QQ plots')
    fig.tight_layout(pad=5.0)
    plt_count = 0
    for phenotype, data in pvals.items():
        ax = axs[plt_count]
        #ax = axs[plt_count // 2, plt_count % 2]
        data[data == 0] = min(data[data != 0])
        data.sort()
        trans_data = -np.log10(data)
        trans_data[trans_data > HEIGHT_CUTOFF] = HEIGHT_CUTOFF
        ax.set_title(f'{phenotype}')
        ax.set_ylabel('-log_10(pval)')
        ax.set_xlabel(f'-log_10(Uniform dist) (total STR loci = {nresults})')
        if phenotype == 'height':
            plot_height_snp_qq(ax)
        ax.plot(
            trans_uniform,
            trans_uniform,
            color=colors['gray'],
            label='y=x, null p-values'
        )
        """
        ax.plot(
            trans_uniform[[0, -1]],
            [-np.log10(GENOME_WIDE_SIG)] * 2,
            color=colors['red'],
            label="p == 8e-5"
        )
        """
        ax.scatter(
            trans_uniform,
            trans_data,
            color=colors['mediumvioletred'],
            marker='.',
            label="STRs on chr2,3"
        )

        crosses_gws = np.argmin(data < GENOME_WIDE_SIG)
        xy = (trans_uniform[crosses_gws],
              trans_data[crosses_gws])
        n_sig_loci = np.sum(data < GENOME_WIDE_SIG)
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

