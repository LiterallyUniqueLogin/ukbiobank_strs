#!/usr/bin/env python3

import os

import numpy as np
import polars as pl
import pandas as pd

ukb  = os.environ['UKB']
workdir = f'{ukb}/post_finemapping/e_splice_STR_overlap'

fdr_thresh = .1

def fdr_cols(p_vals, info_cols):
    p_vals_arr = p_vals.to_numpy()
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
        for count, p_val in enumerate(sort):
            if p_val > (count + 1)/m*fdr_thresh:
                break
        argsort = np.argsort(p_vals)
        findings = ''
        first = True
        for info_arr in info_arrs:
            info_arr[i][info_arr[i] == None] = 'Missing'
            if not first:
                findings += ':'
            else:
                first = False
            findings += info_arr[i][argsort][:count]
        findingss.append(', '.join(findings))
        pss.append(', '.join(str(x) for x in sort[:count]))
    return (pss, findingss, n_tests)

spliceSTR = pl.scan_csv(
    f'{workdir}/yang_spliceSTRs.tab', sep='\t'
).distinct().groupby(
    'hg19_START'
).agg([pl.col('p_values').list(), pl.col('Tissue').list(), pl.col('gene_name').list(), pl.col('str-exon').str.split_exact('-', 1).struct.field('field_1').list().alias('exon')]).collect()

pss, findingss, n_tests = fdr_cols(spliceSTR['p_values'], [spliceSTR['Tissue'], spliceSTR['gene_name'], spliceSTR['exon']])

new_splice = pl.DataFrame({'hg19_START': spliceSTR['hg19_START'], 'splice_p_vals': pd.Series(pss), 'splice_associations (tissue:gene:exonID)': pd.Series(findingss), 'splice_n_tests': n_tests})

eSTR = pl.scan_csv(
    f'{workdir}/yang_eSTRs.tab', sep='\t'
).distinct().groupby(
    'hg19_START'
).agg([pl.col('p_values').list(), pl.col('Tissue').list(), pl.col('gene_name').list()]).collect()

pss, findingss, n_tests = fdr_cols(eSTR['p_values'], [eSTR['Tissue'], eSTR['gene_name']])

new_e = pl.DataFrame({'hg19_START': eSTR['hg19_START'], 'expression_p_vals': pss, 'expression_associations (tissue:gene)': findingss, 'expression_n_tests': n_tests})

results = new_splice.join(new_e, how='outer', on='hg19_START')

filtered = pl.read_csv(f'{workdir}/filtered_STRs.tab', sep='\t').filter(
    pl.col('FILTER') != 'PASS'
).with_column(
    (pl.col('CHROM') + '_' + pl.col('POS').cast(str)).alias('chr_pos')
)

joined = pl.read_csv(
    f'{ukb}/export_scripts/results/causal_STR_candidates_for_publication.tab', sep='\t'
).select(('chr' + pl.col('chrom').cast(str) + '_' + pl.col('start_pos').cast(str)).alias('chr_pos')).distinct().join(
    results,
    how='inner',
    left_on='chr_pos',
    right_on='hg19_START'
).join(
    filtered,
    how='left',
    on='chr_pos'
).filter(
    pl.col('FILTER').is_null()
).drop(['START', 'FILTER', 'CHROM_START', 'CHROM', 'POS']).filter(
    (pl.col('splice_p_vals').str.lengths() > 0) | (pl.col('expression_p_vals').str.lengths() > 0)
).with_columns([
    pl.when(
        pl.col('splice_p_vals').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('splice_p_vals')
    ).alias('splice_p_vals'),
    pl.when(
        pl.col('splice_p_vals').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('splice_associations (tissue:gene:exonID)')
    ).alias('splice_associations (tissue:gene:exonID)'),
    pl.when(
        pl.col('splice_p_vals').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('splice_n_tests')
    ).alias('splice_n_tests'),
    pl.when(
        pl.col('expression_p_vals').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('expression_p_vals')
    ).alias('expression_p_vals'),
    pl.when(
        pl.col('expression_p_vals').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('expression_associations (tissue:gene)')
    ).alias('expression_associations (tissue:gene)'),
    pl.when(
        pl.col('expression_p_vals').str.lengths() == 0
    ).then(
        None
    ).otherwise(
        pl.col('expression_n_tests')
    ).alias('expression_n_tests'),
]).rename({'chr_pos': 'chr_hg19pos'})

joined.write_csv(f'{workdir}/blessed_eSTRs_spliceSTRs.tsv', sep='\t')

