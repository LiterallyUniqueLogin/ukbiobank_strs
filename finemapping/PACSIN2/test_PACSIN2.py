#!/usr/bin/env python3

import glob
import os

import numpy as np

import load_and_filter_genotypes as lfg
import my_regional_gwas

ukb = os.environ['UKB']

varnames = [
    'PACSIN2_compound_repeat_len',
    'PACSIN2_A_STR_partial_len',
    'PACSIN2_TA_STR_len',
    'PACSIN2_CA_STR_len',
    'PACSIN2_T(A|G)_STR_partial_len',
    'PACSIN2_SNP_43385897_C_G',
    'PACSIN2_SNP_43385885_T_C',
    'PACSIN2_SNP_43385893_C_A_SNP_43385924_G_C',
    'PACSIN2_SNP_43385903_C_A',
    'PACSIN2_SNP_43385917_T_G'
]
reflens = [52,2,18,14,18]

def get_gt_itr(samples):
    yield ('motif', 'period', 'ref_len', 'subset_total_per_allele_dosages')
    for varname, reflen in zip(varnames[:5], reflens):
        d = {}
        fnames = glob.glob(f'{ukb}/finemapping/PACSIN2/gts/{varname}*')
        for fname in fnames:
            len_ = int(fname.split('.')[0].split('_')[-1])
            d[len_] = np.load(fname)[samples, :]
        motif = varname.split('_')[1]
        if '|' not in motif:
            period = len(motif)
        else:
            period = 2
        subset_total_per_allele_dosages = {}
        for len_, gts in d.items():
            subset_total_per_allele_dosages[len_] = np.sum(gts)
        yield (d, np.array(sorted(d.keys())), 22, varname, None, (motif, str(period), str(reflen), lfg.dict_str(subset_total_per_allele_dosages)))

for ethnicity in ['south_asian', 'white_brits', 'white_other', 'irish', 'chinese', 'black']:
#for ethnicity in ['white_brits']:
    for phenotype in ['platelet_count', 'mean_platelet_volume', 'platelet_distribution_width']:
            outerdir = f'{ukb}/association/spot_test/{ethnicity}/{phenotype}'
            os.makedirs(outerdir, exist_ok=True)
            #with open('/dev/null', 'w') as outfile:
            with open(f'{outerdir}/PACSIN2.tab', 'w') as outfile:
                my_regional_gwas.perform_regional_gwas_helper(
                    outfile,
                    f'{ukb}/traits/subset_transformed_phenotypes/{ethnicity}/{phenotype}.npy',
                    f'{ukb}/traits/shared_covars/shared_covars.npy',
                    f'{ukb}/traits/phenotypes/{ethnicity}/{phenotype}.npy',
                    get_gt_itr,
                    phenotype,
                    False,
                    None,
                    'strs'
                )
