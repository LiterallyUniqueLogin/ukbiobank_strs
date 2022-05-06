#!/usr/bin/env python3

import ast
import json
import os

import polars as pl

import phenotypes

ukb = os.environ['UKB']

other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']

phenotypes.phenotypes_in_use = ['platelet_distribution_width', 'platelet_count']

def dosages_to_frequencies(dosage_dict_str):
    dosages = ast.literal_eval(dosage_dict_str)
    # drop zero alleles
    dosages = { k: v for k, v in dosages.items() if v != 0 }
    total_dosage = sum(dosages.values())
    return json.dumps({ k: f'{v/total_dosage*100:.2f}%' for k,v in dosages.items() }).replace('"', '')

finemapping_dfs = []
for phenotype in phenotypes.phenotypes_in_use:
    df = pl.scan_csv(
        f'{ukb}/post_finemapping/intermediate_results/finemapping_all_concordance_{phenotype}.tab',
        sep='\t',
        dtypes={
            **{f'{ethnicity}_p_val': float for ethnicity in other_ethnicities},
            **{f'{ethnicity}_coeff': float for ethnicity in other_ethnicities},
            **{f'{ethnicity}_se': float for ethnicity in other_ethnicities}
        }
    ).filter('is_STR')
    fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    with open(fname) as tsv:
        header = tsv.readline().strip()
    assoc_df = pl.scan_csv(
        fname,
        sep='\t',
        skip_rows=1,
        has_header=False,
        with_column_names = lambda _: header.replace('0.05_significance_CI', 'foo', 1).replace('5e-8_significance_CI', 'bar', 1).split('\t') # these duplicate column names won't be used anyway
    ).select([
        'chrom',
        'pos',
        pl.col('subset_total_per_allele_dosages').alias('white_brit_allele_dosages')
    ])
    df = df.join(
        assoc_df,
        how='left',
        on=['chrom', 'pos']
    )
    for ethnicity in other_ethnicities:
        fname = f'{ukb}/association/results_finemapped_only/{ethnicity}/{phenotype}/my_str/results.tab'
        with open(fname) as tsv:
            header = tsv.readline().strip()
        assoc_df = pl.scan_csv(
            fname,
            sep='\t',
            skip_rows=1,
            has_header=False,
            with_column_names = lambda _: header.replace('0.05_significance_CI', 'foo', 1).replace('5e-8_significance_CI', 'bar', 1).split('\t') # these duplicate column names won't be used anyway
        ).select([
            'chrom',
            'pos',
            pl.col('subset_total_per_allele_dosages').alias(f'{ethnicity}_allele_dosages')
        ])
        df = df.join(
            assoc_df,
            how='left',
            on=['chrom', 'pos']
        )
    finemapping_dfs.append(df.collect())
finemapping_results = pl.concat(finemapping_dfs).rename({'pos': 'snpstr_pos'})

finemapping_results = finemapping_results.filter(
    (pl.col('p_val') <= 5e-8) &
    (
        ((pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 0)) | (pl.col('finemap_pip') >= 0.8)
    ).any().over(['chrom', 'snpstr_pos'])
)

pos_table = pl.read_csv(
    f'{ukb}/snpstr/flank_trimmed_vcf/vars.tab',
    sep='\t'
)

finemapping_results = finemapping_results.join(
    pos_table,
    how='left',
    on=['chrom', 'snpstr_pos']
)

repeat_units = pl.read_csv(
    f'{ukb}/snpstr/repeat_units.tab',
    sep='\t'
)

finemapping_results = finemapping_results.join(
    repeat_units,
    how='left',
    on=['chrom', 'snpstr_pos']
)

concordance_cols = pl.read_csv(
    f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_white_blood_cell_count.tab',
    sep='\t',
    n_rows=1
).columns
concordance_results = pl.concat([
    pl.read_csv(
        f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab',
        sep='\t',
        dtypes={col: (float if 'cs' not in col else int) for col in concordance_cols if 'finemap' in col or 'susie' in col or 'p_val' in col}
    ) for phenotype in phenotypes.phenotypes_in_use
    if not os.path.exists(f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab.empty')
]).filter('is_STR').with_column(pl.col('pos').alias('snpstr_pos'))

finemapping_results = finemapping_results.join(
    concordance_results,
    how='left',
    on=['phenotype', 'chrom', 'snpstr_pos']
).with_columns([
    pl.when(pl.col('susie_alpha').is_null()).then(None).when(pl.col('susie_cs') >= 0).then(pl.col('susie_alpha')).otherwise(0).alias('susie_CP'),
    pl.when(pl.col('susie_alpha_hardcall').is_null()).then(None).when(pl.col('susie_cs_hardcall') >= 0).then(pl.col('susie_alpha_hardcall')).otherwise(0).alias('susie_CP_best_guess_genotypes'),
    pl.when(pl.col('susie_alpha_ratio').is_null()).then(None).when(pl.col('susie_cs_ratio') >= 0).then(pl.col('susie_alpha_ratio')).otherwise(0).alias('susie_CP_prior_snps_over_strs'),
    pl.col('finemap_pip').alias('finemap_CP'),
    pl.col('finemap_pip_p_thresh').alias('finemap_CP_pval_thresh_5e-4'),
    pl.col('finemap_pip_mac').alias('finemap_CP_mac_thresh_100'),
    pl.col('finemap_pip_prior_std_derived').alias('finemap_CP_prior_effect_size_0.05%'),
    pl.col('finemap_pip_total_prob').alias('finemap_CP_prior_4_signals'),
    pl.col('finemap_pip_conv_tol').alias('finemap_CP_stopping_thresh_1e-4'),
    pl.col('finemap_pip_ratio').alias('finemap_CP_prior_snps_over_strs'),
    pl.col('finemap_pip_prior_std_low').alias('finemap_CP_prior_effect_size_0.0025%'),
])

finemapping_results = finemapping_results.select([
    'phenotype',
    'chrom',
    pl.col('pos').alias('start_pos'),
    'end_pos',
    pl.col('region').alias('finemapping_region'),
    pl.col('unit').alias('repeat_unit'),
    pl.col('white_brit_allele_dosages').apply(dosages_to_frequencies).alias('white_brit_allele_frequencies'),
    pl.col('p_val').alias('association_p_value'),
    pl.when(pl.col('coeff') > 0).then('+').otherwise('-').alias('direction_of_association'),
    pl.when(
        (pl.col('susie_CP') >= 0.8) &
        (pl.col('finemap_CP') >= 0.8) &
        (pl.col('susie_CP_best_guess_genotypes') >= 0.8) &
        (pl.col('finemap_CP_pval_thresh_5e-4') >= 0.8) &
        (pl.col('finemap_CP_mac_thresh_100') >= 0.8) &
        (pl.col('finemap_CP_prior_effect_size_0.05%') >= 0.8) &
        (pl.col('finemap_CP_prior_4_signals') >= 0.8) &
        (pl.col('finemap_CP_stopping_thresh_1e-4') >= 0.8)
    ).then(
        'confidently'
    ).when(
        (pl.col('susie_CP') >= 0.8) &
        (pl.col('finemap_CP') >= 0.8)
    ).then(
        'doubly'
    ).when(
        (pl.col('susie_CP') >= 0.8) |
        (pl.col('finemap_CP') >= 0.8)
    ).then(
        'singly'
    ).otherwise(
        'not'
    ).alias('finemapping'),
    'susie_CP',
    'finemap_CP',
    'susie_CP_best_guess_genotypes',
    'finemap_CP_pval_thresh_5e-4',
    'finemap_CP_mac_thresh_100',
    'finemap_CP_prior_effect_size_0.05%',
    'finemap_CP_prior_4_signals',
    'finemap_CP_stopping_thresh_1e-4',
    'susie_CP_prior_snps_over_strs',
    'finemap_CP_prior_snps_over_strs',
    'finemap_CP_prior_effect_size_0.0025%',
    pl.sum([
        pl.col(f'{ethnicity}_p_val').cast(str) + pl.lit(', ') for ethnicity in other_ethnicities
    ]).str.replace(', $', '').alias('other_ethnicity_association_p_values'),
    pl.sum([
        pl.when(pl.col(f'{ethnicity}_p_val') > .05).then('NA').when(pl.col('coeff') > 0).then('+').otherwise('-') + pl.lit(', ')
        for ethnicity in other_ethnicities
    ]).str.replace(', $', '').alias('other_ethnicity_effect_directions'),
    *[pl.col(f'{ethnicity}_allele_dosages').apply(dosages_to_frequencies).alias(f'{ethnicity}_allele_frequencies') for ethnicity in other_ethnicities],
])

finemapping_results.write_csv(f'{ukb}/post_finemapping/results/singly_finemapped_strs_for_paper.tab', sep='\t')
finemapping_results.sort(['chrom', 'start_pos']).write_csv(f'{ukb}/post_finemapping/results/singly_finemapped_strs_sorted.tab', sep='\t')

confident_results = finemapping_results.filter(
    (pl.col('finemapping') == 'confidently').any().over(['chrom', 'start_pos']) &
    (pl.col('association_p_val') <= 1e-10)
)
confident_results.write_csv(f'{ukb}/post_finemapping/results/confidently_finemapped_strs_for_paper.tab', sep='\t')
confident_results.sort(['chrom', 'start_pos']).write_csv(f'{ukb}/post_finemapping/results/confidently_finemapped_strs_sorted.tab', sep='\t')
