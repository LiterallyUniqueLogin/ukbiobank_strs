#!/bin/env python3

"""
A simple for loop.

Reads the files using cyvcf2
and pumps them in to scipy
"""

import argparse
from Typing import Tuple
import os

import cyvcf2
import numpy as np
from numpy.typing import ArrayLike
import scipy

ukb = os.environ['UKB']

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


def load_covars():
    """
    Load the sex and population PC covariates.

    Returns
    -------
    Tuple[ArrayLike[int], ArrayLike[int], ArrayLike[float, float]] :
        ids, sex, pcs. The id of each entry in the sex array and each row
        in the pcs array is the corresponding entry in the ids array.

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
    ids = ids_and_sex[:, 0] #each entry here is the id of the corresponding row of pcs
    sex = ids_and_sex[:, 1] #1 = M, 2 = F
    return ids, sex, pcs


def

if __name__ == "__main__":
    main()
