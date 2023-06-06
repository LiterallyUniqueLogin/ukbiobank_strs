#!/usr/bin/env python3

import argparse
import datetime

import numpy as np
import polars as pl

def choose_vars(readme_fname, outcols_fname, str_associations_fname, snp_associations_fname, filter_set_fname, phenotype, chrom, start_pos, end_pos, p_cutoff, mac):
    if mac:
        mac_threshold = int(mac[0])
        snp_mac_fname = mac[1]
        str_mac_fname = mac[2]
        snps_exclude_mac = pl.scan_csv(
            snp_mac_fname,
            sep='\t'
        ).filter(
            pl.col('ALT_CTS') < mac_threshold
        ).select([
            '#POS',
            'REF',
            'ALT'
        ]).collect().pipe(
            lambda df: list(zip(*df.to_dict().values()))
        )

        # need to make that look like a list of strings to polars b/c buggy, so add a single nonsense to it
        snps_exclude_mac.append('asdf')

        strs_exclude_mac = pl.scan_csv(
            str_mac_fname,
            sep='\t'
        ).filter(
            pl.col('mac') < mac_threshold
        ).select(
            'pos'
        ).collect()['pos'].to_list()

    today = datetime.datetime.now().strftime("%Y_%M_%D")
    with open(readme_fname, 'w') as readme:
        readme.write(
            f'Run date: {today}\n'
            f'Choosing variants for which association tests were not skipped and with p <= {p_cutoff}. '
            'SNPs in the filter set are also skipped. '
            f'(Filter set at {filter_set_fname})\n'
        )

    # first choose STRs and SNPs only with p <= p_cutoff to lessen memory burden
    print('Choosing which strs and snps to include ... ', flush=True)
    strs_to_include = set()
    snps_to_include = set()

    if str_associations_fname != "None":
        strs_to_include = pl.scan_csv(
            str_associations_fname,
            sep='\t'
        ).filter(
            (pl.col('chrom') == chrom) &
            (pl.col('pos') >= start_pos) &
            (pl.col('pos') <= end_pos) &
            (pl.col(f'p_{phenotype}') <= p_cutoff)
        ).filter(
            # STRs with duplicate loci that shouldn't have been in the reference panel
            ((pl.col('chrom') != 17) | (pl.col('pos') != 80520458)) &
            ((pl.col('chrom') != 1) | (pl.col('pos') != 247747217)) &
            ((pl.col('chrom') != 1) | (pl.col('pos') != 247848392)) &
            ((pl.col('chrom') != 21) | (pl.col('pos') != 47741815)) &
            ((pl.col('chrom') != 8) | (pl.col('pos') != 145231731))
        ).select(pl.col('pos')).collect().to_numpy().flatten()

        if mac:
            strs_to_include = strs_to_include[~np.isin(strs_to_include, strs_exclude_mac)]

    snps_to_filter = set()
    snps_to_filter = pl.scan_csv(
        filter_set_fname,
        sep='\t'
    ).select([
        pl.col('snp_pos'),
        pl.col('snp_ref'),
        pl.col('snp_alt'),
        pl.lit('1').alias('join_marker')
    ])

    snps_to_include = pl.scan_csv(
        snp_associations_fname,
        sep='\t',
        null_values='NA'
    ).filter(
        (pl.col('#CHROM') == chrom) &
        (pl.col('POS') >= start_pos) &
        (pl.col('POS') <= end_pos) &
        (pl.col('P') <= p_cutoff)
    ).join(
        snps_to_filter,
        how = 'left',
        left_on = ['POS', 'REF', 'ALT'],
        right_on = ['snp_pos', 'snp_ref', 'snp_alt']
    ).filter(
        pl.col('join_marker').is_null()
    ).select([
        'POS',
        'REF',
        'ALT'
    ]).collect().pipe(
        lambda df: list(zip(*df.to_dict().values()))
    )
    # returns a list of tuples

    if mac:
        snps_to_include = [snps_to_include[idx] for idx in np.where(~np.isin(snps_to_include, snps_exclude_mac))[0]]

    snp_sort_tuples = set((pos, 'SNP', ref, alt) for (pos, ref, alt) in snps_to_include)

    if str_associations_fname != "None":
        str_sort_tuples = set((pos, 'STR') for pos in strs_to_include)
        vars_ = snp_sort_tuples.union(str_sort_tuples)

    else:
        vars_ = snp_sort_tuples

    sorted_vars = sorted(vars_)
    sorted_var_names =  [
        f'STR_{tuple[0]}' if tuple[1] == 'STR'
        else f'SNP_{tuple[0]}_{tuple[2]}_{tuple[3]}' if tuple[1] == 'SNP'
        else None # break the sort
        for tuple in sorted_vars
    ]
    assert len(set(sorted_var_names)) == len(sorted_var_names) # make sure is unique

    if str_associations_fname != "None":
        print(f'# STRs: {len(strs_to_include)} # SNPs: {len(snps_to_include)}', flush=True)
    else:
        print(f'Not including STRs. # SNPs: {len(snps_to_include)}', flush=True)

    with open(outcols_fname, 'w') as colfile:
        for var_name in sorted_var_names:
            colfile.write(var_name + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('readme')
    parser.add_argument('outcols')
    parser.add_argument('str_associations')
    parser.add_argument('snp_associations')
    parser.add_argument('filter_set_fname')
    parser.add_argument('phenotype')
    parser.add_argument('chrom', type=int)
    parser.add_argument('start', type=int)
    parser.add_argument('end', type=int)
    parser.add_argument('--threshold', default=5e-4, type=float)
    parser.add_argument('--mac', nargs=3, default=None)
    args = parser.parse_args()

    choose_vars(args.readme, args.outcols, args.str_associations, args.snp_associations, args.filter_set_fname, args.phenotype, args.chrom, args.start, args.end, args.threshold, args.mac)
