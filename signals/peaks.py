#!/usr/bin/env python3

import argparse
import os
import sortedcontainers

import numpy as np
import pandas as pd

import python_array_utils as utils

ukb = os.environ['UKB']

def get_snp_loci(plink_imputed_snp_fname):
    cols = ['#CHROM', 'POS', 'REF', 'ALT', 'P', 'ERRCODE']
    dtypes = utils.get_dtypes(plink_imputed_snp_fname)
    dtypes = {col: dtypes[col] for col in cols}
    csv = utils.df_to_recarray(pd.read_csv(
        plink_imputed_snp_fname,
        header=0,
        delimiter='\t',
        usecols=cols,
        dtype=dtypes
    ))

    # remove xMHC
    csv = csv[~(
        (csv['#CHROM'] == 6) & (25e6 <= csv['POS']) & (csv['POS'] <= 33.5e6)
    )]

    # round to zero past this point
    csv['P'][csv['P'] <= 1e-300] = 0

    csv = csv[csv['ERRCODE'] != 'CONST_OMITTED_ALLELE']

    assert np.all(csv['ERRCODE'] == '.')

    return sortedcontainers.SortedSet(
        iterable = ((row['P'], row['#CHROM'], row['POS'], 'SNP', row['REF'], row['ALT']) for row in iter(csv))
    )



'''
def get_snp_loci(sorted_set, plink_imputed_snp_fname, thresh):
    with open(plink_imputed_snp_fname) as csv:
        next(csv) # skip header
        line_count = 0
        for line in csv:
            line_count += 1
            if line_count % 10000 == 0:
                print(f'{line_count} lines', end='\r', flush=True)
            split = line[:-1].split('\t')
            chrom, pos, ref, alt, p, errcode = (
                split[idx] for idx in (0, 1, 3, 4, 13, 14)
            )
            chrom = int(chrom)
            pos = int(pos)
            if p == 'NA':
                continue
            p = float(p)
       
            # must pass threshold
            if p > thresh:
                continue

            # remove xMHC
            if chrom == 6 and 25e6 <= pos <= 33.5e6:
                continue

            if errcode == 'CONST_OMITTED_ALLELE':
                continue

            # round to zero past this point
            if p <= 1e-300:
                p = 0

            if errcode != '.':
                print(line)
                assert False

            sorted_set.add(
                (p, chrom, pos, 'SNP', ref, alt)
            )
'''

def get_str_loci(phenotype, my_str_fname):
    p_col = f'p_{phenotype}'
    csv = utils.df_to_recarray(pd.read_csv(
        my_str_fname,
        header=0,
        delimiter='\t',
        usecols=['chrom', 'pos', p_col],
        dtype=utils.get_dtypes(my_str_fname)
    ))

    # remove xMHC
    csv = csv[~(
        (csv['chrom'] == 6) & (25e6 <= csv['pos']) & (csv['pos'] <= 33.5e6)
    )]

    # round to zero past this point
    csv[p_col][csv[p_col] <= 1e-300] = 0

    return sortedcontainers.SortedSet(
        iterable = ((row[p_col], row['chrom'], row['pos'], 'STR') for row in iter(csv)),
    )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('spacing', type=int)
    parser.add_argument('my_str_fname')
    parser.add_argument('plink_imputed_snp_fname')
    parser.add_argument('readme_fname')
    parser.add_argument('out_fname')
    args = parser.parse_args()

    phenotype = args.phenotype
    spacing = args.spacing

    print("Loading snps ... ", end='', flush=True)
    snp_loci = get_snp_loci(args.plink_imputed_snp_fname)
    print('done', flush=True)
    print("Loading strs ... ", end='', flush=True)
    str_loci = get_str_loci(phenotype, args.my_str_fname)
    print('done', flush=True)

    # p_val, chrom, pos, type, other
    print("Merging loci ... ", end = '', flush=True)
    loci = snp_loci.union(str_loci)
    print('done', flush=True)

    # chrom, pos, type, p_val, other -> tagged
    peaks = sortedcontainers.SortedDict()

    with open(args.readme_fname, 'w') as readme:
        readme.write(
            f'Considering all variants from input files {args.my_str_fname} and {args.plink_imputed_snp_fname}\n'
            'Choosing association peak variants in the following order (per chromosome): '
            'Round all variants with p < 1e-300 to p=0. '
            'Exclude the extended MHC locus (chr6:25e6-33.5e6). '
            f"Take the variant with the lowest p-value that isn't within {spacing} bp from "
            'any variant already selected, with the following tiebreakers (tiebreakers '
            'should only frequently occur for varaints with p rounded down to 0 and it is '
            'unlikely the second or later tiebreakers will ever be used.):\n'
            '* choose the variant with the smallest starting bp\n'
            '* choose SNPs over STRs\n'
            '* choose SNPs with shorter reference alleles\n'
            '* choose SNPs with lexicographically earlier alternate alleles\n'
            f'This is continued until all variants are examined. '
            f'STR/SNP peak varaints within {spacing} of variants of the other type that pass the '
            'p-value threshold are marked as being tagged by the other type of variants.'
        )

    print('Gathering peaks ... ', flush=True)
    for locus in loci:
        potential_peak = (*locus[1:4], locus[0], *locus[4:])
        idx = peaks.bisect_left(potential_peak)
        too_close = False
        if idx > 0:
            lower_peak, _ = peaks.peekitem(idx-1)
            if locus[1:3] == ('chr1', '25250919'):
                print(lower_peak)
            if lower_peak[0] == potential_peak[0] and lower_peak[1] >= potential_peak[1] - spacing:
                too_close = True
                if lower_peak[2] != potential_peak[2]:
                    peaks[lower_peak] = True
        if idx < len(peaks):
            upper_peak, _ = peaks.peekitem(idx)
            if locus[1:3] == ('chr1', '25250919'):
                print(upper_peak)
            if upper_peak[0] == potential_peak[0] and upper_peak[1] <= potential_peak[1] + spacing:
                too_close = True
                if upper_peak[2] != potential_peak[2]:
                    peaks[upper_peak] = True
        if too_close:
            continue
        peaks[potential_peak] = False
        print(f'Found a peak! {len(peaks)} total peaks', end='\r', flush=True)
    print('\n')
    print('Done gathering peaks', flush=True)

    print('Writing out peaks ... ', flush=True, end='')
    with open(args.out_fname, 'w') as outfile:
        outfile.write('chrom\tpos\tvariant_type\tp_value\ttagged_by_other_variant_type\tref_(snp_only)\talt_(snp_only)\n')
        outfile.flush()
        for peak, tagged in peaks.items():
            str_peak = [str(item) for item in peak]
            if len(str_peak) == 5:
                str_peak.append('')
                str_peak.append('')
            outfile.write('\t'.join(str_peak[:4]) + '\t' + str(tagged) + '\t' + '\t'.join(str_peak[4:]) + '\n')
            outfile.flush()
    print('done', flush=True)

if __name__ == '__main__':
    main()
