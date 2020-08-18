import argparse
import os
from matplotlib.colors import CSS4_COLORS as colors
import matplotlib.pyplot as plt
import numpy as np

ukb = os.environ['UKB']

def plot_entropy_by_locus(ax):
    chr_lens = np.loadtxt(
        f'{ukb}/misc_data/genome/chr_lens.txt',
        usecols=[1],
        skiprows=1
    )
    prior_len = 0
    ax.set_title('Locus vs entropy')
    ax.set_xlabel('Genomic position, in bp, lining up each chr')
    ax.set_ylabel('Entropy')
    for chrom in range(1, 23):
        if chrom > 1:
            prior_len += chr_lens[chrom-2]
        if chrom % 2 == 0:
            color = colors['mediumvioletred']
        else:
            color = colors['deeppink']
        data = np.loadtxt(
            f'{ukb}/side_analyses/entropy/full_genome/chr{chrom}.tab',
            usecols=[1, 3],
            skiprows=1
        )
        ax.scatter(data[:, 0] + prior_len, data[:, 1], color=color, marker='.')

def main():
    fig, ax = plt.subplots(1, figsize=(30, 5))
    plot_entropy_by_locus(ax)
    plt.savefig(f'{ukb}/side_analyses/entropy/full_genome/entropy_by_locus.png')


if __name__ == "__main__":
    main()
