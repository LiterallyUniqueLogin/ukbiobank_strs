#!/usr/bin/env python3

import argparse

import numpy as np
import pandas as pd

# default spacing is 250kb between variants (125kb at each side of the leftmost and right most variants)
# which is same spacing that plink uses in its clumping algorithm (though its algorithm
# is more nuanced)
def generate_clumps(chr_lens_fname, iters, spacing=250000):
    chr_lens = np.genfromtxt(
        chr_lens_fname,
        usecols=[1],
        skip_header=1,
        dtype=int
    )
    curr_items = [next(itr) for itr in iters]
    # currently only designed for integer chromosome numbers

    curr_chrom, curr_pos = min(curr_items)
    curr_pos_min = max(curr_pos - spacing//2, 0)
    for idx, val in enumerate(curr_items):
        if val == (curr_chrom, curr_pos):
            curr_items[idx] = next(iters[idx], (np.inf, np.inf))
    while True:
        next_chrom, next_pos = min(curr_items)
        if next_chrom < curr_chrom:
            raise ValueError("Chroms out of order")
        if next_chrom == np.inf:
            if curr_pos_min != np.inf:
                yield curr_chrom, curr_pos_min, min(curr_pos + spacing//2, chr_lens[int(curr_chrom) - 1])
            break

        if next_chrom != curr_chrom or curr_pos + spacing < next_pos:
            yield curr_chrom, curr_pos_min, min(curr_pos + spacing//2, chr_lens[int(curr_chrom) - 1])
            curr_chrom = next_chrom
            curr_pos_min = max(next_pos - spacing//2, 0)

        curr_pos = next_pos

        for idx, val in enumerate(curr_items):
            if val == (next_chrom, next_pos):
                curr_items[idx] = next(iters[idx], (np.inf, np.inf))

def plink_snp_output_itr(plink_snp_fname, thresh):
    csv = pd.read_csv(
        plink_snp_fname,
        header=0,
        delimiter='\t',
        usecols=['#CHROM', 'POS', 'P']
    ).to_numpy()

    csv = csv[csv[:, 2]  < thresh, :]

    for linenum in range(csv.shape[0]):
        yield tuple(csv[linenum, :2])

def my_str_output_itr(phenotype, my_str_fname, thresh):
    csv = pd.read_csv(
        my_str_fname,
        header=0,
        delimiter='\t',
        usecols=['chrom', 'pos', f'p_{phenotype}']
    ).to_numpy()

    csv = csv[csv[:, 2]  < thresh, :]

    for linenum in range(csv.shape[0]):
        yield tuple(csv[linenum, :2])

def filter_MHC_itr(itr):
    for chrom, pos in itr:
        if chrom == 6 and 25e6 <= pos <= 33.5e6:
            continue
        yield (chrom, pos)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('chr_lens_fname')
    parser.add_argument('plink_imputed_snp_fname')
    parser.add_argument('out_fname')
    parser.add_argument('--remove-skips', action='store_true')
    parser.add_argument('--my-str-fname')
    args = parser.parse_args()

    phenotype = args.phenotype

    pthresh = 5e-8
    snp_itr = plink_snp_output_itr(args.plink_imputed_snp_fname, pthresh)

    if args.my_str_fname:
        results = pd.read_csv(
            args.my_str_fname,
            header=0,
            usecols=['chrom', 'pos', f'p_{phenotype}', 'locus_filtered'],
            delimiter='\t',
            encoding='UTF-8',
            dtype={'chrom': int, 'pos': int, f'p_{phenotype}': float, 'locus_filtered': str}
        )
        results = results.loc[results[f'p_{phenotype}'] <= 5e-4, :]
        results = results.loc[results['locus_filtered'] == 'False', :]

        itrs = [
            snp_itr,
            my_str_output_itr(phenotype, args.my_str_fname, pthresh)
        ]
    else:
        itrs = [snp_itr]

    itrs = [filter_MHC_itr(itr) for itr in itrs]
    with open(args.out_fname[:-4] + '_README.txt', 'w') as readme:
        readme.write(
            f'Clumping results from my_str and plink_snp runs for phenotype {phenotype} '
            'by cenetering a 250kb interval (125kb in each direction) around each variant '
            f'from either run that passes the {pthresh:g} threshold and then joining all '
            'overlapping intervals.\n'
        )
    with open(args.out_fname, 'w') as outfile:
        outfile.write('chrom\tstart\tend\tany_strs\n')
        outfile.flush()
        prev_clump = None
        for clump in generate_clumps(args.chr_lens_fname, itrs):
            clump = [int(val) for val in clump]
            if prev_clump is not None and prev_clump[0] == clump[0] and prev_clump[2] >= clump[1]:
                print(prev_clump, clump)
                assert False
            if args.my_str_fname:
                strs = np.any(
                    (results['chrom'] == clump[0]) &
                    (clump[1] <= results['pos']) &
                    (results['pos'] <= clump[2])
                )
            else:
                strs = False
            if args.remove_skips and ((
                phenotype == 'urate' and
                clump[0] == 4 and
                clump[1] == 8165642 and
                clump[2] == 11717761
            ) or (
                phenotype == 'total_bilirubin' and
                clump[0] == 12 and
                clump[1] == 19976272 and
                clump[2] == 22524428
            ) or (
                phenotype == 'alkaline_phosphatase' and
                clump[0] == 1 and
                clump[1] == 19430673 and
                clump[2] == 24309348
            )):
                continue
            outfile.write(f'{clump[0]}\t{clump[1]}\t{clump[2]}\t{strs}\n')
            outfile.flush()
            prev_clump = clump

if __name__ == '__main__':
    main()
