#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd
import sortedcontainers

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

    # round to zero past this point
    csv['P'][csv['P'] <= 1e-300] = 0

    csv = csv[csv['ERRCODE'] != 'CONST_OMITTED_ALLELE']

    assert np.all(csv['ERRCODE'] == '.')

    return sortedcontainers.SortedSet(
        iterable = ((row['P'], row['#CHROM'], row['POS'], 'SNP', row['REF'], row['ALT']) for row in iter(csv))
    )

def get_str_loci(phenotype, my_str_fname):
    p_col = f'p_{phenotype}'
    csv = utils.df_to_recarray(pd.read_csv(
        my_str_fname,
        header=0,
        delimiter='\t',
        usecols=['chrom', 'pos', p_col],
        dtype=utils.get_dtypes(my_str_fname)
    ))

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

    print('Gathering potential peaks ... ', flush=True)
    for locus in loci:
        potential_peak = (*locus[1:4], locus[0], *locus[4:])
        idx = peaks.bisect_left(potential_peak)
        if idx > 0:
            lower_peak, _ = peaks.peekitem(idx-1)
            if lower_peak[0] == potential_peak[0] and lower_peak[1] >= potential_peak[1] - spacing:
                continue
        if idx < len(peaks):
            upper_peak, _ = peaks.peekitem(idx)
            if upper_peak[0] == potential_peak[0] and upper_peak[1] <= potential_peak[1] + spacing:
                continue
        # add a new peak, mark as untagged for now
        peaks[potential_peak] = np.nan

    print('Removing slope elements, tagging peaks ... ', flush=True)
    for locus in loci:
        locus_p_val = locus[0]
        locus_as_peak = (*locus[1:4], locus[0], *locus[4:])
        idx = peaks.bisect_left(locus_as_peak)
        del_lower = False
        if idx > 0:
            lower_peak, _ = peaks.peekitem(idx-1)
            lower_peak_p_val = lower_peak[3]
            if (
                # check for nearness
                lower_peak[0] == locus_as_peak[0] and
                lower_peak[1] >= locus_as_peak[1] - spacing
            ):
                if lower_peak_p_val > locus_p_val:
                    print('del lower')
                    # Remove slope element
                    del_lower = True
                elif lower_peak[2] != locus_as_peak[2]:
                    # tag actual peaks
                    if np.isnan(peaks[lower_peak]):
                        peaks[lower_peak] = locus_p_val
                    else:
                        peaks[lower_peak] = min(peaks[lower_peak], locus_p_val)

        del_upper = False
        if idx < len(peaks):
            upper_peak, _ = peaks.peekitem(idx)
            upper_peak_p_val = upper_peak[3]
            if (
                # check for nearness
                upper_peak[0] == locus_as_peak[0] and
                upper_peak[1] - spacing <= locus_as_peak[1]
            ):
                if upper_peak_p_val > locus_p_val:
                    print('del upper')
                    # Remove slope element
                    del_upper = True
                elif upper_peak[2] != locus_as_peak[2]:
                    # tag actual peaks
                    if np.isnan(peaks[upper_peak]):
                        peaks[upper_peak] = locus_p_val
                    else:
                        peaks[upper_peak] = min(peaks[upper_peak], locus_p_val)

        if del_lower:
            del peaks[lower_peak]
        if del_upper:
            del peaks[upper_peak]


    for peak, tag_p_val in peaks.items():
        if not np.isnan(tag_p_val) and peak[3] > tag_p_val:
            print(peak, tag_p_val)
            assert False

    print('Done gathering peaks', flush=True)
    print('Writing out peaks ... ', flush=True, end='')
    with open(args.out_fname, 'w') as outfile:
        outfile.write('chrom\tpos\tvariant_type\tp_value\tp_value_other_variant_type\tref_(snp_only)\talt_(snp_only)\n')
        outfile.flush()
        for peak, tagged in peaks.items():
            if peak[1] == 1 and peak[2] == 9341128:
                print('printing peak')
            str_peak = [str(item) for item in peak]
            if len(str_peak) == 5:
                str_peak.append('')
                str_peak.append('')
            outfile.write('\t'.join(str_peak[:4]) + '\t' + str(tagged) + '\t' + '\t'.join(str_peak[4:]) + '\n')
            outfile.flush()
    print('done', flush=True)

if __name__ == '__main__':
    main()
