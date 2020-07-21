#!/bin/env python3

"""
A simple for loop.

Reads the files using cyvcf2
and pumps them in to scipy
"""

import argparse
import os

import cyvcf2 
import numpy as np
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

    pcs = np.loadtxt(
        f'{ukb}/misc_data/EGA/ukb_sqc_v2.txt'
        usecols=range(-43, -3)
    ) #156MB
    ids_and_sex = np.loadtxt(
        f'{ukb}/microarray/ukb46122_cal_chr1_v2_s488282.fam',
        usecols=(0, 4),
        dtype=int,
    ) 
    ids = ids_and_sex[:,0] #each entry here is the id of the corresponding row of pcs
    sex = ids_and_sex[:,1] #1 = M, 2 = F
    

if __name__ == "__main__":
    main()
