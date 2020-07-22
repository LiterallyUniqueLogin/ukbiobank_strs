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
    pcs = np.loadtxt(
        f'{ukb}/misc_data/EGA/ukb_sqc_v2.txt',
        usecols=range(-43, -3)
    ) #156MB

    ids_and_sex = np.loadtxt(
        f'{ukb}/microarray/ukb46122_cal_chr1_v2_s488282.fam',
        usecols=(0, 4),
        dtype=int,
    )
    np_df = np.concatenate((ids_and_sex, pcs), axis=1)
    return pd.DataFrame(
        np_df,
        columns=("id", "sex", *(f"pc{col}" for col in range(1, 41)))
    ).astype({"id" : "int", "sex" : "int"}).astype({"id": "category"})

    # for sex 1 = M, 2 = F


def load_height(covar_df):
    """
    Append a height column to the covariate df.

    Height is measured in cm.

    Returns
    -------
    pd.DataFrame :
       Same as covar_df, but with two rows appended.
       The first is 'height' as a float
       (heights are measured in half cms)
       The second is 'height_sampling'
       and is categorical 0, 1 or 2
       specifying which visit height was retrieved on.
       Only first visit where height was retrieved is
       recorded.
    """
    height = np.loadtxt(
        f'{ukb}/main_dataset/extracted_data/height.txt',
        delimiter=" ",
        skiprows=4
    )
    height_df = pd.DataFrame(
        height,
        columns = ["id", "height", "height_sampling"]
    ).astype({"id" : "int", "height_sampling" : "category"}).astype({"id": "category"})
    return pd.merge(
        covar_df,
        height_df,
        how="left",
        on="id"
    )


def match_ids(pheno_ids, pheno, covar_ids, covars):
    """
    Return a subset of the datasets with the shared IDs.

    The dataset will be in the covar_ids order.

    Parameters
    ----------
    pheno_ids :
        A list of ids
    pheno :
        An array where each id corresponds to a row
    covar_ids :
        A list of ids
    covars :
        A list of arrays where, for each array, each id
        corresponds to a row

    Returns
    -------
    Tuple[ArrayLike[int], ArrayLike, ArrayLike :
        The three arrays the shared ids, the phenotypes for those ids and the
        covars for those ids
    """

    # reduce covar to shared ids
    covar_id_comp = np.equal(
        pheno_ids.reshape(1, -1),
        covar_ids.reshape(-1, 1)
    )
    covar_subset = np.any(covar_id_comp, axis=1)
    covar_ids_shared = covar_ids[covar_subset]
    covars_shared = [covar[covar_subset, :] for covar in covars]

    # reduce pheno to shared ids and sort
    # to covars id order
    pheno_ids_comp = np.equal(
        pheno_ids.reshape(1, -1),
        covar_ids_shared.reshape(-1, 1)
    )
    pheno_subset = np.any(pheno_ids_comp, axis=0)
    pheno_sort = np.argmax(pheno_ids_comp, axis=0)
    # SnO standards for shared and ordered
    pheno_ids_SnO = pheno_ids[pheno_subset][pheno_sort[pheno_subset]]
    pheno_SnO = pheno[pheno_subset, :][pheno_sort[pheno_subset], :]
    assert np.all(np.equals(pheno_ids_SnO, covar_ids_shared))
    return (
        covar_ids_shared,
        pheno_SnO,
        covars_shared
    )


def main():  # noqa: D103
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('filtering_run_name')
    parser.add_argument('imputation_run_name')
    parser.add_argument('samples_file')
    args = parser.parse_args()

    filtering_run_name = args.filtering_run_name
    imputation_run_name = args.imputation_run_name
    samples_file = args.samples_file

	vcf_floc = (f'{ukb}/str_imputed/runs/{imputation_run_name}/'
				f'vcfs/strs_only/chr{chrom}.vcf.gz')
    """

    df = load_covars()
    print(df.head())
    return
    unmatched_height = load_height()
    print(unmatched_height[:5, :3])

    ids, height, covars = match_ids(
        unmatched_height[:, 0],
        unmatched_height[:, 1],
        unmatched_covar_ids,
        [unmatched_sex.reshape(-1, 1), unmatched_pcs]
    )
    print(ids[:5], height[:5, :3], [covar[:5, :3] for covar in covars])


if __name__ == "__main__":
    main()

