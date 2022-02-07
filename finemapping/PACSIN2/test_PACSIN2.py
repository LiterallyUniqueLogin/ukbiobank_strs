#!/usr/bin/env python3

import os

import my_regional_gwas
import load_PACSIN2

ukb = os.environ['UKB']

if __name__ == '__main__':
    for ethnicity in ['south_asian', 'white_brits', 'white_other', 'irish', 'chinese', 'black']:
        for phenotype in ['platelet_count', 'mean_platelet_volume', 'platelet_distribution_width']:
                outerdir = f'{ukb}/association/spot_test/{ethnicity}/{phenotype}'
                os.makedirs(outerdir, exist_ok=True)
                with open(f'{outerdir}/PACSIN2.tab', 'w') as outfile:
                    my_regional_gwas.perform_regional_gwas_helper(
                        outfile,
                        f'{ukb}/traits/subset_transformed_phenotypes/{ethnicity}/{phenotype}.npy',
                        f'{ukb}/traits/shared_covars/shared_covars.npy',
                        f'{ukb}/traits/phenotypes/{ethnicity}/{phenotype}.npy',
                        load_PACSIN2.get_gt_itr,
                        phenotype,
                        False,
                        None,
                        'strs'
                    )
