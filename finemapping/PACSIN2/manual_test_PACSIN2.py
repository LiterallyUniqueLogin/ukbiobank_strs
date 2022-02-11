#!/usr/bin/env python3

import os

import numpy as np
from statsmodels.regression.linear_model import OLS

import load_PACSIN2
import python_array_utils
import sample_utils

ukb = os.environ['UKB']


scovs = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')
print('phenotype\tethnicity\tcompound coeff single test (s.d.)\tA coeff (s.d.)\tTA coeff (s.d.)\tCA coeff (s.d.)')
for phenotype in ('mean_platelet_volume', 'platelet_distribution_width', 'platelet_count'):
    for ethnicity in ('white_brits', 'black', 'south_asian', 'chinese', 'irish', 'white_other'):
        print(f'{phenotype}\t{ethnicity}', end='')
        pcovs = np.load(f'{ukb}/traits/subset_transformed_phenotypes/{ethnicity}/{phenotype}.npy')
        samps = sample_utils.get_ordered_samples_phenotype(ethnicity, phenotype).reshape(-1, 1)
        covs = python_array_utils.merge_arrays(python_array_utils.merge_arrays(samps, pcovs), scovs)

        outcomes = covs[:, 1]
        covs = covs[:, 2:]

        samp_idx = sample_utils.get_samples_idx_phenotype(ethnicity, phenotype)
        itr = load_PACSIN2.get_gt_itr(samp_idx)
        next(itr) # skip details
        compound, a, ta, ca = [next(itr)[0] for _ in range(4)]
        compound, a, ta, ca = [np.sum([len_*np.sum(dosages, axis=1) for len_, dosages in dosage_dict.items()], axis=0).reshape(-1, 1) for dosage_dict in (compound, a, ta, ca)]
        stds = [np.std(gts) for gts in (compound, a, ta, ca)]
        compound, a, ta, ca = [(gts - np.mean(gts))/np.std(gts) for gts in (compound, a, ta, ca)]

        for gt in a, ta, ca:
            result = OLS(
                outcomes,
                np.hstack((covs, np.ones((covs.shape[0], 1)), gt))
            ).fit()
            print(f'\t{result.pvalues[-1]}', end='')
            #print(f'\t{result.params[-1]/stds[0]:.2g} ({result.bse[-1]/stds[0]:.2g})', end='')
        '''

        result = OLS(
            outcomes,
            np.hstack((covs, np.ones((covs.shape[0], 1)), a, ta, ca))
        ).fit()
        for idx in (-3, -2, -1):
            print(f'\t{result.params[idx]/stds[idx]:.2g} ({result.bse[idx]/stds[idx]:.2g})', end='')
        print(flush=True)
        '''
        
