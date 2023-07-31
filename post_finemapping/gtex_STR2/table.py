#!/usr/bin/env python3

import argparse

import numpy as np
import polars as pl
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('confidently_finemapped_table')
parser.add_argument('QTL_results')
parser.add_argument('--methylation', action='store_true', default=False)

args = parser.parse_args()

bf_thresh = .05

def bf_cols(p_vals, info_cols):
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
            if p_val <= bf_thresh/m:
                success = True
                thresh = m - count
                break
        if not success:
            thresh = 0
        argsort = np.argsort(p_vals)
        findings = ''
        first = True
        if 'cg04111102' in list(info_cols[3][i].to_numpy()):
            print(info_cols[3][i].to_numpy())
            print('----')
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

trait_assocs = pl.read_csv(
    args.confidently_finemapped_table,
    sep='\t'
).filter(
    pl.col('finemapping') == 'confidently'
).select([
    ('chr' + pl.col('chrom').cast(str) + '_' + pl.col('start_pos (hg19)').cast(str)).alias('chrom_pos19'),
    ('chr' + pl.col('chrom').cast(str) + '_' + pl.col('start_pos (hg38)').cast(str)).alias('chrom_pos38'),
    *[('chr' + pl.col('chrom').cast(str) + '_' + (pl.col('start_pos (hg38)') + offset).cast(str)).alias(f'chrom_pos38_{offset}') # prepare to do a hacky 10-tolerance asof join
      for offset in range(-10, 11)],
    'phenotype',
    pl.col('association_p_value').cast(str)
]).groupby('chrom_pos19').agg([
    pl.col('chrom_pos38').first(),
    *[pl.col(f'chrom_pos38_{offset}').first() for offset in range(-10, 11)],
    pl.col('association_p_value').list(),
    pl.col('phenotype').list()
]).with_columns([
    pl.col('association_p_value').arr.join(','),
    pl.col('phenotype').arr.join(',')
])

if not args.methylation:
    col_name = 'str-gene'
else:
    col_name = 'STR_DNAmID'

qtl_str = pl.scan_csv(
    args.QTL_results,
    sep='\t'
).with_columns([
    pl.col(col_name).str.split_exact('-', 1).struct.field('field_0').alias('hg38'),
    pl.when(
        pl.col('slope') > 0
    ).then(
        pl.lit('+')
    ).otherwise(
        pl.lit('-')
    ).alias('effect_direction')
])

if args.methylation:
    qtl_str = qtl_str.with_columns([
        pl.lit('WholeBlood').alias('Tissue'),
        pl.col(col_name).str.split_exact('-', 1).struct.field('field_1').alias('gene_name')
    ])

# shrink memory usage by dropping unnecessary columns
qtl_str = qtl_str.select([
    col_name,
    'p_values',
    'effect_direction',
    'Tissue',
    'gene_name',
    'hg38'
]).collect()

print('got here', flush=True)

qtl_str = pl.concat([
    qtl_str.join(
        trait_assocs,
        left_on='hg38',
        right_on=f'chrom_pos38_{offset}'
    ).drop([
        f'chrom_pos38_{offset2}' for offset2 in range(-10, 11) if offset2 != offset
    ])
    for offset in range(-10, 11)
])

print('concat succeeded', flush=True)
print(qtl_str.shape, qtl_str.columns, flush=True)

#qtl_str = qtl_str.unique().groupby(
qtl_str = qtl_str.groupby(
    'chrom_pos19'
).agg([
    pl.col('phenotype').first(),
    pl.col('association_p_value').first(),
    pl.col('p_values').list(),
    pl.col('effect_direction').list(),
    pl.col('Tissue').list(),
    pl.col('gene_name').list(),
    pl.col(col_name).str.split_exact('-', 1).struct.field('field_1').list().alias('target')
])

print('grouping succeeded', flush=True)
print(qtl_str.shape, qtl_str.columns, flush=True)

print(qtl_str.filter(pl.col('target').arr.contains('cg04111102')))
exit()

pss, findingss, n_tests = bf_cols(qtl_str['p_values'], [qtl_str['effect_direction'], qtl_str['Tissue'], qtl_str['gene_name'], qtl_str['target']])

#print(f'----- eSTR ---------')
#for str_, ps, findings in zip(qtl_str['chrom_pos19'], pss, findingss):
#    assert len(ps.split(',')) == len(findings.split(','))
#    if len(ps) == 0:
#        continue
#    finding_list = [x.strip() for x in findings.split(',')]
#    #if fname == 'eSTR':
#    finding_list = [x[:(len(x) - x[::-1].index(':') - 1)] for x in finding_list]
#    finding_list = [x[2:] + ':' + x[0] for x in finding_list]
#    print(str_, ', '.join(f'{finding} ({p})' for p, finding in zip((x.strip() for x in ps.split(',')), finding_list)))

total_qtl_str = pl.DataFrame({
    'chrom_pos': qtl_str['chrom_pos19'],
    'phenotype': qtl_str['phenotype'],
    'phenotype_association_p_value': qtl_str['association_p_value'],
    'expression_association_p_vals': pd.Series(pss),
    'expression_associations (tissue:target)': pd.Series(findingss),
    'n_tests': n_tests
}).sort('chrom_pos')

total_qtl_str = total_qtl_str.filter(
    pl.col('expression_association_p_vals').str.lengths() > 0
).sort('chrom_pos').select(['chrom_pos', pl.all().exclude('^chrom_pos$')])

if not args.methylation:
    out = 'qtl_STRs.tab'
else:
    out = 'meQTL_STRs.tab'
total_qtl_str.to_pandas().to_csv(out, sep='\t', index=False)
