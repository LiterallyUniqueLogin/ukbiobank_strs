#!/usr/bin/env python3

import argparse
import datetime
import os

import polars as pl

ukb = os.environ['UKB']

def load_gts(readme_fname, outcols_fname, phenotype, chrom, start_pos, end_pos, use_PACSIN2):
    if use_PACSIN2:
        assert int(chrom) == 22

    filter_set_fname = f'{ukb}/finemapping/str_imp_snp_overlaps/chr{chrom}_to_filter.tab'

    p_cutoff = 5e-4

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

    strs_to_include = pl.scan_csv(
        f'{ukb}/association/results/{phenotype}/my_str/results.tab',
        sep='\t'
    ).filter(
        (pl.col('chrom') == chrom) &
        (pl.col('pos') >= start_pos) &
        (pl.col('pos') <= end_pos) &
        (pl.col(f'p_{phenotype}') <= p_cutoff)
    ).select(pl.col('pos')).collect().to_numpy().flatten()

    assert len(strs_to_include) != 0

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
        f'{ukb}/association/results/{phenotype}/plink_snp/results.tab',
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
        pl.col('POS'),
        pl.col('REF'),
        pl.col('ALT')
    ]).collect().pipe(
        lambda df: list(zip(*df.to_dict().values()))
    )
    # returns a list of tuples

    snp_sort_tuples = set((pos, 'SNP', ref, alt) for (pos, ref, alt) in snps_to_include)
    str_sort_tuples = set((pos, 'STR') for pos in strs_to_include)
    vars_ = snp_sort_tuples.union(str_sort_tuples)
    if use_PACSIN2:
        vars_.remove((43385872, 'STR'))
        vars_.add((43385866, 'PACSIN2_STR'))
        vars_.add((43385875, 'PACSIN2_STR'))
        vars_.add((43385893, 'PACSIN2_STR'))
    sorted_vars = sorted(vars_)
    sorted_var_names =  [
        f'STR_{tuple[0]}' if tuple[1] == 'STR'
        else f'SNP_{tuple[0]}_{tuple[2]}_{tuple[3]}' if tuple[1] == 'SNP'
        else f'PACSIN2_STR_{tuple[0]}' if tuple[1] == 'PACSIN2_STR'
        else None # break the sort
        for tuple in sorted_vars
    ]
    assert len(set(sorted_var_names)) == len(sorted_var_names) # make sure is unique

    print(f'# STRs: {len(strs_to_include)} # SNPs: {len(snps_to_include)}', flush=True)

    with open(outcols_fname, 'w') as colfile:
        for var_name in sorted_var_names:
            colfile.write(var_name + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('readme')
    parser.add_argument('outcols')
    parser.add_argument('phenotype')
    parser.add_argument('chrom')
    parser.add_argument('start')
    parser.add_argument('end')
    parser.add_argument('--three-PACSIN2-STRs', action='store_true', default=False)
    args = parser.parse_args()

    load_gts(args.readme, args.outcols, args.phenotype, args.chrom, args.start, args.end, args.three_PACSIN2_STRs)
