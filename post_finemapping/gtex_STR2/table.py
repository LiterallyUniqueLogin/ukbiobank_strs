#!/usr/bin/env python3

import os

import numpy as np
import polars as pl
import pandas as pd

ukb = os.environ['UKB']

fdr_thresh = .05/3

def fdr_cols(p_vals, info_cols):
    p_vals_arr = p_vals.to_pandas().to_numpy()
    info_arrs = []
    for info_col in info_cols:
        info_arrs.append(info_col.to_numpy())
    pss = []
    n_tests = []
    findingss = []
    for i in range(len(p_vals_arr)):
        ps = []
        p_vals = p_vals_arr[i]
        sort = np.sort(p_vals)
        m = len(p_vals)
        n_tests.append(m)
        success = False
        for count, p_val in enumerate(sort[::-1]):
            if p_val <= (m - count)/m*fdr_thresh:
                success = True
                thresh = m - count
                break
        if not success:
            thresh = 0
        argsort = np.argsort(p_vals)
        findings = ''
        first = True
        for info_arr in info_arrs:
            info_arr[i][info_arr[i] == None] = 'Missing'
            if not first:
                findings += ':'
            else:
                first = False
            findings += info_arr[i][argsort][:thresh]
        findingss.append(', '.join(findings))
        pss.append(', '.join(f'{x:.3g}' for x in sort[:thresh]))
    return (pss, findingss, n_tests)

coords = []
for fname in  'snpstr_strs_19.bed', 'snpstr_strs_38.bed':
    coords.append(pl.read_csv(
        f'{ukb}/side_analyses/exome_strs/intermediate_files/{fname}',
        sep='\t',
        header=None,
        new_columns=['chrom', 'start', 'end', '_', '__', 'str']
    ).with_columns([
        (pl.col('chrom') + '_' + (pl.col('start') + 1 + offset).cast(str)).alias(f'chrom_pos_{offset}')
        for offset in range(-10, 11)
    ]))
coords_df = coords[0].join(
    coords[1], on=['chrom', 'str'], suffix='_38'
).select([
    'chrom_pos_0',
    *[f'chrom_pos_{offset}_38' for offset in range(-10, 11)]
]).rename({'chrom_pos_0': 'chrom_pos'})

trait_assocs = pl.read_csv(
    f'{ukb}/post_finemapping/results/confidently_finemapped_strs_for_paper.tab',
    sep='\t'
).filter(
    pl.col('finemapping') == 'confidently'
).select([
    ('chr' + pl.col('chrom').cast(str) + '_' + pl.col('start_pos').cast(str)).alias('chrom_pos'),
    'phenotype',
    pl.col('association_p_value').cast(str)
]).groupby('chrom_pos').agg([
    pl.col('association_p_value').list(),
    pl.col('phenotype').list()
]).with_columns([
    pl.col('association_p_value').arr.join(','),
    pl.col('phenotype').arr.join(',')
])

coords_df = coords_df.join(
    trait_assocs,
    left_on='chrom_pos',
    right_on='chrom_pos'
)
#print(coords_df.shape)

qtl_strs = []
yang_dir = '/expanse/projects/gymreklab/yal084_storage/share_with_Jonathan'
for fname, col_name in ('eSTR', 'str-gene'), ('STR', 'str-exon'), ('eISOFORM', 'str-isoform'):
    qtl_str = pl.read_csv(
        f'{yang_dir}/{fname}_GB_650pc_combined_fdr10p.csv',
        sep='\t'
    ).with_column(
        pl.col(col_name).str.split_exact('-', 1).struct.field('field_0').alias('hg38')
    )
    qtl_str = pl.concat([
        qtl_str.join(
            coords_df,
            left_on='hg38',
            right_on=f'chrom_pos_{offset}_38'
        ).drop([
            f'chrom_pos_{offset2}_38' for offset2 in range(-10, 11) if offset2 != offset
        ])
        for offset in range(-10, 11)
    ]).with_column(
        pl.when(
            pl.col('slope') > 0
        ).then(
            pl.lit('+')
        ).otherwise(
            pl.lit('-')
        ).alias('effect_direction')
    )
    print(qtl_str.shape)
    exit()

    qtl_str = qtl_str.distinct().groupby(
        'chrom_pos'
    ).agg([pl.col('phenotype').first(), pl.col('association_p_value').first(), pl.col('p_values').list(), pl.col('effect_direction').list(), pl.col('Tissue').list(), pl.col('gene_name').list(), pl.col(col_name).str.split_exact('-', 1).struct.field('field_1').list().alias('target')])
    #print(qtl_str.shape)

    pss, findingss, n_tests = fdr_cols(qtl_str['p_values'], [qtl_str['effect_direction'], qtl_str['Tissue'], qtl_str['gene_name'], qtl_str['target']])

    print(f'----- {fname} ---------')
    for str_, ps, findings in zip(qtl_str['chrom_pos'], pss, findingss):
        assert len(ps.split(',')) == len(findings.split(','))
        if len(ps) == 0:
            continue
        finding_list = [x.strip() for x in findings.split(',')]
        if fname == 'eSTR':
            finding_list = [x[:(len(x) - x[::-1].index(':') - 1)] for x in finding_list]
        finding_list = [x[2:] + ':' + x[0] for x in finding_list]
        print(str_, ', '.join(f'{finding} ({p})' for p, finding in zip((x.strip() for x in ps.split(',')), finding_list)))

    qtl_strs.append(pl.DataFrame({'chrom_pos': qtl_str['chrom_pos'],'phenotype':qtl_str['phenotype'], 'association_p_value':qtl_str['association_p_value'],
                                  'p_vals': pd.Series(pss), 'associations (tissue:target)': pd.Series(findingss), 'n_tests': n_tests}).sort('chrom_pos'))

total_qtl_str = qtl_strs[0].join(
    qtl_strs[1].drop(['phenotype', 'association_p_value']),
    how='outer',
    on='chrom_pos',
    suffix='_splice'
).join(
    qtl_strs[2].drop(['phenotype', 'association_p_value']),
    how='outer',
    on='chrom_pos',
    suffix='_isoform'
).rename({'p_vals': 'p_vals_expression', 'associations (tissue:target)' : 'associations (tissue:target)_expression', 'n_tests' : 'n_tests_expression'})

total_qtl_str = total_qtl_str.filter(
    (pl.col('p_vals_expression').str.lengths() > 0) | (pl.col('p_vals_splice').str.lengths() > 0) | (pl.col('p_vals_isoform').str.lengths() > 0)
).with_columns([
    pl.when(
        pl.col('p_vals_expression').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('p_vals_expression')
    ).alias('p_vals_expression'),
    pl.when(
        pl.col('p_vals_expression').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('associations (tissue:target)_expression')
    ).alias('associations (tissue:target)_expression'),
    pl.when(
        pl.col('p_vals_expression').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('n_tests_expression')
    ).alias('n_tests_expression'),
    pl.when(
        pl.col('p_vals_splice').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('p_vals_splice')
    ).alias('p_vals_splice'),
    pl.when(
        pl.col('p_vals_splice').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('associations (tissue:target)_splice')
    ).alias('associations (tissue:target)_splice'),
    pl.when(
        pl.col('p_vals_splice').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('n_tests_splice')
    ).alias('n_tests_splice'),
    pl.when(
        pl.col('p_vals_isoform').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('p_vals_isoform')
    ).alias('p_vals_isoform'),
    pl.when(
        pl.col('p_vals_isoform').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('associations (tissue:target)_isoform')
    ).alias('associations (tissue:target)_isoform'),
    pl.when(
        pl.col('p_vals_isoform').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('n_tests_isoform')
    ).alias('n_tests_isoform'),
]).sort('chrom_pos').select(['chrom_pos', pl.all().exclude('^chrom_pos$')])

total_qtl_str.to_pandas().to_csv('blessed_qtl_STRs.tab', sep='\t', index=False)
