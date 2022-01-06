#!/usr/bin/env python3

import argparse
import json
import os

import numpy as np
import polars as pl

min_ld_thresh = 0.9

def load_dir(pheno, region, dir_):
    with open(f'{dir_}/converged.txt') as converged:
        assert converged.read().strip() == 'TRUE'

    alphas = pl.scan_csv(
        f'{dir_}/alpha.tab',
        sep='\t',
        has_header=False
    ).collect().to_numpy().T
    susie_pips=1-np.prod(1-alphas, axis=1)

    df = pl.scan_csv(
        f'{dir_}/colnames.txt',
        has_header=False,
        with_column_names=lambda _: ['var_name']
    ).with_column(
        pl.lit(1).alias('row_number')
    ).with_columns([
        pl.col('row_number').cumsum(),
        pl.lit(None, int).alias('cs_num'),
        pl.lit(region).alias('region'),
        pl.lit(pheno).alias('phenotype'),
        pl.Series(susie_pips).alias('susie_pip'),
        pl.lit(None, float).alias('susie_cs_pip')
    ])

    for cs_num in range(50):
        cs_num+=1
        cs_fname = f'{dir_}/cs{cs_num}.txt'
        if not os.path.exists(cs_fname):
            continue
        with open(cs_fname) as cs:
            var_nums = [int(var_num) for var_num in next(cs).strip().split()]
            next(cs)
            min_ld = float(next(cs).split()[0])
            if min_ld < min_ld_thresh:
                continue
            df = df.with_columns([
                pl.when(
                    pl.col('row_number').is_in(var_nums)
                ).then(
                    pl.when(
                        ~pl.col('cs_num').is_null()
                    ).then(-1).otherwise(
                        cs_num
                    )
                ).otherwise(
                    pl.col('cs_num')
                ).alias('cs_num'),
                pl.when(
                    pl.col('row_number').is_in(var_nums)
                ).then(
                    pl.Series(alphas[:, cs_num-1])
                ).otherwise(
                    pl.col('susie_cs_pip')
                ).alias('susie_cs_pip')
            ])

    df = df.with_column(
        pl.when(
            pl.col('cs_num') != -1
        ).then(
            pl.col('susie_cs_pip')
        ).otherwise(
            -1
        ).alias('susie_cs_pip')
    )
    df = df.filter(
        pl.col('var_name').str.contains('^STR') &
        ~pl.col('cs_num').is_null() &
        (pl.col('susie_pip') > 0.05)
    ).drop('row_number')

    return df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('json_dirs_list_fname')
    parser.add_argument('outfname')
    args = parser.parse_args()

    with open(args.json_dirs_list_fname) as json_dirs_list_file:
        dirs_list = json.loads(json_dirs_list_file.read().strip())

    dfs = []
    for pheno, region, dir_ in dirs_list:
        if pheno == 'urate' and region == '4_8165642_11717761':
            continue
        print(f'Loading region {pheno}: {region}', flush=True)
        dfs.append(load_dir(pheno, region, dir_))
    pl.concat(dfs).with_columns([
        pl.col('susie_pip').round(4),
        pl.col('susie_cs_pip').round(4),
    ]).collect().to_csv(
        args.outfname, sep='\t'
    )

if __name__ == '__main__':
    main()
