#!/usr/bin/env python3

import argparse
import os

import cyvcf2
import numpy as np
from statsmodels.regression.linear_model import OLS

import python_array_utils
import sample_utils

ukb = os.environ['UKB']

# do an association test of the combined dosage of the listed alleles vs the listed phenotypes in each ethnicity
parser = argparse.ArgumentParser()
parser.add_argument('chrom', type=int)
parser.add_argument('pos', type=int)
parser.add_argument('--phenotypes', nargs='+')
parser.add_argument('--alleles', type=int, nargs='+') #allele indicies, i.e. ths SNP is present in alleles 0 (ref), 2 (2nd alt) and 5
args = parser.parse_args()

scovs = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')
var = next(cyvcf2.VCF(f'{ukb}/str_imputed/runs/first_pass/vcfs/annotated_strs/chr{args.chrom}.vcf.gz')(f'{args.chrom}:{args.pos}'))
print('phenotype\tethnicity\tp-val')
for phenotype in args.phenotypes:
    for ethnicity in ('white_brits', 'black', 'south_asian', 'chinese', 'irish', 'white_other'):
        print(f'{phenotype}\t{ethnicity}', end='')
        pcovs = np.load(f'{ukb}/traits/subset_transformed_phenotypes/{ethnicity}/{phenotype}.npy')
        samps = sample_utils.get_ordered_samples_phenotype(ethnicity, phenotype).reshape(-1, 1)
        covs = python_array_utils.merge_arrays(python_array_utils.merge_arrays(samps, pcovs), scovs)

        outcomes = covs[:, 1]
        covs = covs[:, 2:]

        samp_idx = sample_utils.get_samples_idx_phenotype(ethnicity, phenotype)
        gts = np.zeros(covs.shape[0])
        for allele in args.alleles:
            if allele > 0:
                gts += (var.format('AP1') + var.format('AP2'))[samp_idx, allele]
            else:
                assert allele == 0
                gts += 1 - np.sum((var.format('AP1') + var.format('AP2'))[samp_idx, :], axis=1)

        if np.all(gts == 0):
            print('\t1')
            continue
        gts = (gts - np.mean(gts))/np.std(gts)

        result = OLS(
            outcomes,
            np.hstack((covs, np.ones((covs.shape[0], 1)), gts.reshape(-1, 1)))
        ).fit()
        print(f'\t{result.pvalues[-1]}')

