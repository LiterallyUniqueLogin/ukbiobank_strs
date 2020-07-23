#!/bin/env python3

"""
A simple for loop.

Reads the files using cyvcf2
and pumps them in to scipy
"""

import argparse
import os

import cyvcf2
import pandas as pd
import numpy as np
import scipy

ukb = os.environ['UKB']


def load_covars():
    """
    Load the sex and population PC covariates.

    Returns
    -------
    pd.DataFrame :
        A categorical column 'id' which serves as the logical index
        An integer column 'sex' with M = 1 and F = 2
        Float columns pc1 ... pc40 for the population structure pcs

    Notes
    -----
    If interested in the source of this data, read the READMEs in the
    directories of the loaded files.
    """
    print("Loading covariates ... ", end="", flush=True)
    pcs = pd.read_csv(
        f'{ukb}/misc_data/EGA/ukb_sqc_v2.txt',
        sep=" ",
        header=None,
        names=list(f"pc{col}" for col in range(1, 41)),
        usecols=range(25, 65),
        dtype=float
    ) #156MB

    ids_and_sex = pd.read_csv(
        f'{ukb}/microarray/ukb46122_cal_chr1_v2_s488282.fam',
        usecols=(0, 4),
        header=None,
        names=['id', 'sex'],
        dtype=int,
        sep=" "
    )
    try:
        return pd.concat((ids_and_sex, pcs), axis=1)
    finally:
        print("done")


def load_height(df):
    """
    Append height columns to the df.

    Height is measured in cm.

    Returns
    -------
    pd.DataFrame :
       Same as the input, but with two rows appended.
       The first is 'height' as a float
       (heights are measured in half cms)
       The second is 'height_sampling'
       and is categorical 0, 1 or 2
       specifying which visit height was retrieved on.
       Only first visit where height was retrieved is
       recorded.
    """
    print("Loading height ... ", end="", flush=True)
    height_df = pd.read_csv(
        f'{ukb}/main_dataset/extracted_data/height.txt',
        names=["id", "height", "height_sampling"],
        dtype={"id": "int",
               "height": float,
               "height_sampling": "category"},
        skiprows=4,
        sep=" "
    )

    try:
        return pd.merge(
            df,
            height_df,
            how="left",
            on="id"
        )
    finally:
        print("done")


def load_bilirubin(df):
    """
    Append bilirubin columns to the df.

    Measured in umol/L

    Returns
    -------
    pd.DataFrame :
       Same as input df, but with two float rows appended.
       They are 'direct bilirubin' and
       'total bilirubin'. 'indirect bilirubin'
       can be calculated as total - direct.
       These only reflect measurements made at the first
       visit, follow up measurements are ignored
       as I do not know if bilirubin values change as people
       age or if I should be adjusting for that.
       Unsure why someone would get a value for
       total or direct bilirubin but not the other,
       so for now make sure both exist or NaN them both
    """
    print("Loading bilirubin ... ", end="", flush=True)
    bilirubin_df = pd.read_csv(
        f'{ukb}/main_dataset/extracted_data/bilirubin.csv',
        names=["id", "direct bilirubin", "total bilirubin"],
        header=0,
        usecols=[0, 1, 3],
        dtype={"id": "int",
               "direct bilirubin" : float,
               "total bilirubin" : float},
        quotechar='"'
    )
    direct_nan = np.isnan(bilirubin_df['direct bilirubin'])
    total_nan = np.isnan(bilirubin_df['total bilirubin'])
    bilirubin_df.loc[direct_nan, 'total bilirubin'] = np.nan
    bilirubin_df.loc[total_nan, 'direct bilirubin'] = np.nan

    """
    both_nan = np.logical_and(direct_nan, total_nan)
    print("Direct nan count", np.sum(direct_nan))  # 103892
    print("total nan count", np.sum(total_nan))  # 34945
    print("both nan count", np.sum(both_nan))  # 34601
    """

    try:
        return pd.merge(
            df,
            bilirubin_df,
            how="left",
            on="id"
        )
    finally:
        print("done")


def main():  # noqa: D103
    parser = argparse.ArgumentParser()
    parser.add_argument('filtering_run_name')
    parser.add_argument('imputation_run_name')
    args = parser.parse_args()

    filtering_run_name = args.filtering_run_name
    imputation_run_name = args.imputation_run_name

    """
	vcf_floc = (f'{ukb}/str_imputed/runs/{imputation_run_name}/'
				f'vcfs/strs_only/chr{chrom}.vcf.gz')
    """

    df = load_covars()
    df = load_height(df)
    df = load_bilirubin(df)

    """
    height_nan = np.isnan(df['height'])
    bilirubin_nan = np.isnan(df['direct bilirubin'])
    both_nan = np.logical_and(height_nan, bilirubin_nan)
    print("height nan count", np.sum(height_nan))  # 1487
    print("bilirubin nan count", np.sum(bilirubin_nan))  # 93578
    print("both nan count", np.sum(both_nan))  # 437
    """

    # For ease of use, only work with samples which are
    # present for both
    # Fix later
    df = df[~np.isnan(df['height'])]
    df = df[~np.isnan(df['total bilirubin'])]
    # Already calculated the largest unrelated subset of a filtering list
    # based on height. Will need to do this more intelligently when we
    # don't use the same samples for bilirubin and height
    sample_subset = pd.read_csv(
        f'{ukb}/sample_qc/runs/{filtering_run_name}/combined_unrelated.sample',
        usecols=[0],
        names=['id'],
        dtype=int,
        skiprows=1,
        sep=" "
    )
    df = pd.merge(df, sample_subset, how='inner', on='id')
 
    df = df.astype({"id": "category"})

    # Run the association assuming all phenotypes are present simultaneously

if __name__ == "__main__":
    main()

