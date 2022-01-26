#!/usr/bin/env python3

import argparse
import os
import os.path

import polars as pl

ukb = os.environ['UKB']

def different_vars(finemap_vars_to_compare_fname, phenotype, chrom, start_pos, end_pos):
    filter_set_fname = f'{ukb}/finemapping/str_imp_snp_overlaps/chr{chrom}_to_filter.tab'

    p_cutoff = 5e-2

    # first choose STRs and SNPs only with p <= p_cutoff to lessen memory burden
    #print('Loading strs and snps var list ... ', flush=True)
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
    var_names =  {
        f'STR_{tuple[0]}' if tuple[1] == 'STR' else f'SNP_{tuple[0]}_{tuple[2]}_{tuple[3]}'
        for tuple in vars_
    }

    with open(finemap_vars_to_compare_fname) as finemap_vars_to_compare_file:
        # ignore header
        next(finemap_vars_to_compare_file)
        finemap_vars = {line.split()[0] for line in finemap_vars_to_compare_file.readlines() if line.split()[0]}

    if finemap_vars == var_names:
        return None
    else:
        return list(x for x in var_names if x not in finemap_vars), list(x for x in finemap_vars if x not in var_names)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('chrom')
    parser.add_argument('start_pos')
    parser.add_argument('end_pos')
    args = parser.parse_args()
    phenotype = args.phenotype
    chrom = args.chrom
    start = args.start_pos
    end = args.end_pos

    with open(f'post_finemapping/intermediate_results/finemap_var_discrep_{phenotype}_{chrom}_{start}_{end}.txt', 'w') as output:
        finemap_colnames_file = f'finemapping/finemap_results/{phenotype}/{chrom}_{start}_{end}/finemap_input.z'
        if not os.path.exists(finemap_colnames_file):
            output.write('FINEMAP: missing colnames\n')
            exit()
        finemap_diff = different_vars(finemap_colnames_file, phenotype, chrom, start, end)
        if finemap_diff:
            output.write(f'FINEMAP: {" ".join(finemap_diff[0])} , {" ".join(finemap_diff[1])}\n')
