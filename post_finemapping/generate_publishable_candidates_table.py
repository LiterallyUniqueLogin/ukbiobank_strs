#!/usr/bin/env python3

import argparse

import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outloc')
parser.add_argument('table_before_concordance')
parser.add_argument('concordant_table')
parser.add_argument('start_poses_table')
args = parser.parse_args()

table_before_concordance = pl.read_csv(
    args.table_before_concordance,
    sep='\t'
).with_columns([
    pl.col(f'{ethnicity}_population_allele_frequencies').str.replace_all("'", '')
    for ethnicity in ['black', 'south_asian', 'chinese', 'irish', 'white_other']
]).rename({
    f'{ethnicity}_population_allele_frequencies': f'{ethnicity}_allele_frequencies'
    for ethnicity in ['black', 'south_asian', 'chinese', 'irish', 'white_other']
}).with_column(
    pl.col('relation_to_gene').str.replace('intergenic;', '')
)

concordant_table =  pl.read_csv(
    args.concordant_table,
    sep='\t'
).rename({
    'region': 'finemapping_region'
}).drop([
    'p_val',
    'susie_alpha',
    'finemap_pip',
    'susie_alpha_hardcall',
    'finemap_pip_total_prob',
    'finemap_pip_prior_std_derived',
    'finemap_pip_conv_tol',
    'finemap_pip_mac',
    'finemap_pip_gt_thresh'
]).rename({
    'susie_alpha_ratio': 'susie_pip_ratio_4',
    'finemap_pip_ratio': 'finemap_pip_ratio_4',
    'finemap_pip_prior_std_low': 'finemap_pip_prior_std_0.005'
})

concordant_table = concordant_table.with_columns([
    pl.Series(
        [f'{p:.2f}' for p in concordant_table['finemap_pip_prior_std_0.005'].to_numpy()]
    ).alias('finemap_pip_prior_std_0.005'),
    pl.Series(
        [f'{p:.2f}' for p in concordant_table['finemap_pip_ratio_4'].to_numpy()]
    ).alias('finemap_pip_ratio_4'),
    pl.Series(
        [f'{p:.2f}' for p in concordant_table['susie_pip_ratio_4'].to_numpy()]
    ).alias('susie_pip_ratio_4')
])

start_poses_table = pl.read_csv(
    args.start_poses_table,
    sep='\t'
).rename({'pos': 'actual_pos'})

concordant_table = concordant_table.join(
    start_poses_table,
    how='left',
    left_on=['chrom', 'pos'],
    right_on=['chrom', 'snpstr_pos']
).drop('pos').rename({'actual_pos': 'pos'})

table_before_concordance = table_before_concordance.join(
    concordant_table,
    how='inner',
    left_on=['phenotype', 'chrom', 'start_pos'],
    right_on=['phenotype', 'chrom', 'pos']
).rename({
    'other_ethnic_association_ps': 'other_ethnicity_association_p_values',
    'end_pos': 'end_pos_copy',
    'finemapping_region': 'finemapping_region_copy'
}).drop([
    'FINEMAP_pcausal',
    'SuSiE_CS_pcausal',
    'transcribed'
])

table_before_concordance = table_before_concordance.with_column(pl.Series(
    [f'{p:.2e}' for p in table_before_concordance['association_p_value'].to_numpy()]
).alias('association_p_value'))

table_before_concordance.insert_at_idx(
    1 + table_before_concordance.find_idx_by_name('start_pos'),
    table_before_concordance['end_pos_copy'].alias('end_pos')
)

table_before_concordance.insert_at_idx(
    1 + table_before_concordance.find_idx_by_name('end_pos'),
    table_before_concordance['finemapping_region_copy'].alias('finemapping_region')
)

table_before_concordance.drop([
    'finemapping_region_copy',
    'end_pos_copy',
    'multiallelicness',
    'SuSiE_pcausal'
]).to_csv(args.outloc, sep='\t')
