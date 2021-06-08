#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd

ukb = os.environ['UKB']

chr_lens = np.genfromtxt(
    f'{ukb}/misc_data/genome/chr_lens.txt',
    usecols=[1],
    skip_header=1,
    dtype=int
)

# default spacing is 250kb which is same spacing that
# plink uses in its clumping algorithm (though its algorithm
# is more nuanced)
def generate_clumps(iters, spacing=250000):
    curr_items = [next(itr) for itr in iters]
    # currently only designed for integer chromosome numbers

    curr_chrom, curr_pos = min(curr_items)
    curr_pos_min = max(curr_pos - spacing, 0)
    for idx, val in enumerate(curr_items):
        if val == (curr_chrom, curr_pos):
            curr_items[idx] = next(iters[idx], (np.inf, np.inf))
    while True:
        next_chrom, next_pos = min(curr_items)
        if next_chrom < curr_chrom:
            raise ValueError("Chroms out of order")
        if next_chrom == 6 and 25e6 <= next_pos <= 33.5e6:
            # skip MHC
            continue
        if next_chrom == np.inf:
            if curr_pos_min != np.inf:
                yield curr_chrom, curr_pos_min, min(curr_pos + spacing, chr_lens[int(curr_chrom) - 1])
            break

        if next_chrom != curr_chrom or curr_pos + spacing < next_pos:
            yield curr_chrom, curr_pos_min, min(curr_pos + spacing, chr_lens[int(curr_chrom) - 1])
            curr_chrom = next_chrom
            curr_pos_min = max(next_pos - spacing, 0)

        curr_pos = next_pos

        for idx, val in enumerate(curr_items):
            if val == (next_chrom, next_pos):
                curr_items[idx] = next(iters[idx], (np.inf, np.inf))

def plink_snp_output_itr(phenotype, thresh):
    csv = pd.read_csv(
        f'{ukb}/association/results/{phenotype}/plink_snp/results.tab',
        header=0,
        delimiter='\t',
        usecols=['#CHROM', 'POS', 'P']
    ).to_numpy()

    csv = csv[csv[:, 2]  <= thresh, :]

    for linenum in range(csv.shape[0]):
        yield tuple(csv[linenum, :2])

def my_str_output_itr(phenotype, thresh):
    csv = pd.read_csv(
        f'{ukb}/association/results/{phenotype}/my_str/results.tab',
        header=0,
        delimiter='\t',
        usecols=['chrom', 'pos', f'p_{phenotype}']
    ).to_numpy()

    csv = csv[csv[:, 2]  <= thresh, :]

    for linenum in range(csv.shape[0]):
        yield tuple(csv[linenum, :2])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    args = parser.parse_args()
    phenotype = args.phenotype
    pthresh = 5e-8
    itrs = [plink_snp_output_itr(phenotype, pthresh), my_str_output_itr(phenotype, pthresh)]
    with open(f'{ukb}/finemapping/signal_clumps/{phenotype}_README.txt', 'w') as readme:
        readme.write(
            f'Clumping results from my_str and plink_snp runs for phenotype {phenotype} '
            'by cenetering a 500kb interval (250kb in each direction) around each variant '
            f'from either run that passes the {pthresh:g} threshold and then joining all '
            'overlapping intervals.\n'
        )
    with open(f'{ukb}/finemapping/signal_clumps/{phenotype}.tab', 'w') as outfile:
        outfile.write('chrom\tstart\tend\n')
        outfile.flush()
        for clump in generate_clumps(itrs):
            outfile.write(f'{int(clump[0])}\t{int(clump[1])}\t{int(clump[2])}\n')
            outfile.flush()

if __name__ == '__main__':
    main()
