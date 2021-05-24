#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd

ukb = os.environ['UKB']

def p_filter_iter(itr, thresh):
    '''
    Parameters
    ----------
    itr:
        yields tuples (chrom, pos, p_val)
    '''

    for chrom, pos, p_val in itr:
        if p_val <= thresh:
            yield (chrom, pos)

# default spacing is 250kb which is same spacing that
# plink uses in its clumping algorithm (though its algorithm
# is more nuanced)
def generate_clumps(iters, pthresh = 5e-8, spacing=250000):
    iters = [p_filter_iter(itr, pthresh) for itr in iters]
    curr_items = [next(itr) for itr in iters]
    # currently only designed for integer chromosome numbers

    curr_chrom, curr_pos = min(curr_items)
    curr_pos_min = curr_pos - spacing
    for idx, val in enumerate(curr_items):
        if val == (curr_chrom, curr_pos):
            curr_items[idx] = next(iters[idx], (np.inf, np.inf))
    while True:
        next_chrom, next_pos = min(curr_items)
        if next_chrom == np.inf:
            break
        for idx, val in enumerate(curr_items):
            if val == (next_chrom, next_pos):
                curr_items[idx] = next(iters[idx], (np.inf, np.inf))
        if next_chrom < curr_chrom:
            raise ValueError("Chroms out of order")
        if next_chrom != curr_chrom or curr_pos + spacing < next_pos:
            yield (curr_chrom, curr_pos_min, curr_pos + spacing)
            curr_chrom = next_chrom
            curr_pos_min = next_pos - spacing
            curr_pos = next_pos
        else:
            curr_pos = next_pos

def plink_snp_output_itr(phenotype):
    csv = pd.read_csv(
        f'{ukb}/association/results/{phenotype}/plink_snp/results.tab',
        header=0,
        delimiter='\t',
        usecols=['#CHROM', 'POS', 'P']
    ).to_numpy()

    for linenum in range(csv.shape[0]):
        yield csv[linenum, :]

def my_str_output_itr(phenotype):
    csv = pd.read_csv(
        f'{ukb}/association/results/{phenotype}/my_str/results.tab',
        header=0,
        delimiter='\t',
        usecols=['chrom', 'pos', f'p_{phenotype}']
    ).to_numpy()

    for linenum in range(csv.shape[0]):
        yield csv[linenum, :]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    args = parser.parse_args()
    phenotype = args.phenotype
    itrs = [plink_snp_output_itr(phenotype), my_str_output_itr(phenotype)]
    with open(f'{ukb}/finemapping/signal_clumps/{phenotype}_README.txt', 'w') as readme:
        readme.write(
            f'Clumping results from my_str and plink_snp runs for phenotype {phenotype} '
            'by cenetering a 500kb interval (250kb in each direction) around each variant '
            'from either run that passes the 5e-8 threshold and then joining all '
            'overlapping intervals.\n'
        )
    with open(f'{ukb}/finemapping/signal_clumps/{phenotype}.txt', 'w') as outfile:
        outfile.write('chrom\tstart\tend\n')
        outfile.flush()
        for clump in generate_clumps(itrs):
            outfile.write(f'{int(clump[0])}\t{int(clump[1])}\t{int(clump[2])}\n')
            outfile.flush()

if __name__ == '__main__':
    main()
