#!/usr/bin/env python3

import os

import polars as pl
import numpy as np

ukb = os.environ['UKB']

paths = [
    'old',
    'n_configs_100_rep1',
    'n_configs_100_rep2',
    'n_configs_fixed_rep1',
    'n_configs_fixed_rep2',
    'stable_attempt_rep1',
    'stable_attempt_rep2',
]
paths = [f'{ukb}/temp/finemap_probs/first_pass_df_{path}.tab' for path in paths]

for comp, path_idx1, path_idx2, in [
    ('old v new', 0, 1),
    ('two reps 100', 1, 2),
    ('100 vs fixed', 2, 3),
    ('two reps fixed', 3, 4),
    ('stable attempt vs fixed', 4, 5),
    ('two reps stable attempt', 5, 6)
]:
    print(f'\n---------{comp}----------')
    dfs = []
    for path in path_idx1, path_idx2:
        df = pl.read_csv(
            paths[path],
            sep='\t'
        ).select([
            'chrom', 'pos', 'is_STR', 'varname', 'p_val', 'susie_alpha',
            ((pl.col('susie_alpha') >= 0.8) & (pl.col('p_val') <= 5e-8)).alias('singly_finemapped')
        ])
        assert df.shape == df.unique(subset=['chrom', 'varname']).shape

        dfs.append(df)

    merged = dfs[0].join(
        dfs[1],
        how='inner',
        on=['chrom', 'varname']
    ).with_column(
        (pl.col('susie_alpha') - pl.col('susie_alpha_right')).abs().alias('diff')
    )
    assert abs(merged.shape[0] == dfs[0].shape[0]) <= 5, (merged.shape[0], dfs[0].shape[0])
    assert abs(merged.shape[0] == dfs[1].shape[0]) <= 5, (merged.shape[0], dfs[1].shape[0])

    for exp, cond in (pl.col('is_STR'), 'STRs'), (~pl.col('is_STR'), 'SNPs'):
        print(f'----{cond}:')
        curr_merged = merged.filter(exp)
        print(
            'n shared finemapped: ',
            curr_merged.filter(pl.col('singly_finemapped') & pl.col('singly_finemapped_right')).shape[0]
        )
        print(
            'diffs: ',
            list(np.sort(curr_merged.filter(pl.col('singly_finemapped') & pl.col('singly_finemapped_right'))['diff'].to_numpy())[::-1])
        )
        print(
            'n finemapped only left: ',
            curr_merged.filter(pl.col('singly_finemapped') & ~pl.col('singly_finemapped_right')).shape[0]
        )
        print(
            'diffs: ',
            list(np.sort(curr_merged.filter(pl.col('singly_finemapped') & ~pl.col('singly_finemapped_right'))['diff'].to_numpy())[::-1])
        )
        print(
            'n finemapped only right: ',
            curr_merged.filter(~pl.col('singly_finemapped') & pl.col('singly_finemapped_right')).shape[0]
        )
        print(
            'diffs: ',
            list(np.sort(curr_merged.filter(~pl.col('singly_finemapped') & pl.col('singly_finemapped_right'))['diff'].to_numpy())[::-1]),
            flush=True
        )

