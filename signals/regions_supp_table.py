#!/usr/bin/env python3

import os

import polars as pl

import phenotypes

ukb = os.environ['UKB']

dfs = []
for phenotype in phenotypes.phenotypes_in_use:
    dfs.append(pl.scan_csv(
        f'{ukb}/signals/regions/{phenotype}.tab',
        sep='\t'
    ).with_column(pl.lit(phenotype).alias('phenotype')))

pl.concat(dfs).collect().with_column(
    (
        (
            (pl.col('phenotype') == 'total_bilirubin') &
            (pl.col('chrom') == 12) &
            (pl.col('start') == 19976272) &
            (pl.col('end') == 22524428)
        ) |
        (
            (pl.col('phenotype') == 'urate') &
            (pl.col('chrom') == 4) &
            (pl.col('start') == 8165642) &
            (pl.col('end') == 11717761)
        ) |
        (
            (pl.col('phenotype') == 'alkaline_phosphatase') &
            (pl.col('chrom') == 1) &
            (pl.col('start') == 19430673) &
            (pl.col('end') == 24309348)
        )
    ).alias('filtered_due_to_computation_burden')
).select(
    ['phenotype', 'chrom', 'start', 'end', 'filtered_due_to_computation_burden']
).to_csv(f'{ukb}/export_scripts/results/supp_table_2_finemapping_regions.tab', sep='\t')
