import argparse
import os

from matplotlib.colors import CSS4_COLORS as colors
import matplotlib.pyplot as plt
import numpy as np


ukb = os.environ['UKB']
GENOME_WIDE_SIG = 5e-8


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_name")
    args = parser.parse_args()

    results = np.loadtxt(
        f'{ukb}/association/runs/{args.run_name}/results.txt',
        skiprows=1,
        usecols=(0, 1, 4, 6, 8, 10)
    )
    chroms = results[:, 0]
    pos = results[:, 1]
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
        trans_data[trans_data > 20] = 20
        ax.set_title(f'{phenotype}')
        ax.set_ylabel('-log_10(pval)')
        ax.set_xlabel(f'-log_10(Uniform dist) (total loci = {nresults})')
        ax.plot(
            trans_uniform,
            trans_uniform,
            color=colors['gray'],
            label='y=x, null p-valus'
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
            label="qq plot"
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
