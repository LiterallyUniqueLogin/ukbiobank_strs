#!/usr/bin/env python3

import argparse
import os
import os.path
import time

import numpy as np
import pandas as pd

import python_array_utils as utils

ukb = os.environ['UKB']

def append_mfi_to_plink(result_fname, out_fname):
    # Load plink SNP results
    print(f"Loading plink SNP results at {result_fname}  ... ", end='', flush=True)
   
    start_time = time.time()

    plink_results = pd.read_csv(
        result_fname,
        sep='\t',
        usecols=(0,1,2,3,4,13,14),
        encoding='UTF-8',
        header=0,
        names=('chr', 'pos', 'id', 'ref', 'alt', 'p_val', 'error')
    )

    print("and appending maf/info data ... ", end='', flush=True)
    mfi_dfs = []
    text_squash = 0
    for chrom in np.unique(plink_results['chr']):
        print(text_squash*' ' + text_squash*'\b', end='')
        text = f'(loading chr {chrom})'
        text_squash = len(text)
        print(text + text_squash*'\b', end='', flush=True)
        df = pd.read_csv(
            f'{ukb}/array_imputed/ukb_mfi_chr{chrom}_v3.txt',
            sep='\t',
            usecols=(1,2,3,4,5,7),
            header=0,
            names=('id', 'pos', 'ref', 'alt', 'maf', 'info')
        )
        df['chr'] = chrom
        mfi_dfs.append(df)
    mfi_df = pd.concat(mfi_dfs, axis=0)

    print(text_squash*' ' + text_squash*'\b', end='')
    text = '(merging and converting dataframes)'
    text_squash = len(text)
    print(text + text_squash*'\b', end='', flush=True)
    plink_results = utils.df_to_recarray(plink_results.merge(mfi_df, how='left', on=['chr', 'pos', 'id', 'ref', 'alt']))

    np.save(out_fname, plink_results)
    assert np.all(~np.isnan(plink_results['maf'])) and np.all(~np.isnan(plink_results['info']))

    print(text_squash*' ' + text_squash*'\b', end='')
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("result_fname")
    parser.add_argument("out_fname")
    args = parser.parse_args()

    append_mfi_to_plink(args.result_fname, args.out_fname)

if __name__ == "__main__":
    main()

