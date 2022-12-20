import os

import numpy as np

import python_array_utils as utils

ukb = os.environ['UKB']

def get_all_samples():
    imp_snp_samples_fname = f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample'
    with open(imp_snp_samples_fname) as imp_snp_samples_file:
        imp_snp_samples = np.array([line.split()[0] for line in imp_snp_samples_file][2:], dtype=int).reshape(-1, 1)

    hardcall_snp_samples_fname = f'{ukb}/microarray/ukb46122_hap_chr1_v2_s487314.sample'
    with open(hardcall_snp_samples_fname) as hardcall_snp_samples_file:
        hardcall_snp_samples = np.array([line.split()[0] for line in hardcall_snp_samples_file][2:], dtype=int).reshape(-1, 1)

    # TODO this seems wrong
    compare_idxs = (hardcall_snp_samples > 0) & (imp_snp_samples > 0)
    assert np.all(hardcall_snp_samples[compare_idxs] == imp_snp_samples[compare_idxs])
    return hardcall_snp_samples

def samples_array_with_indicator(sample_fname):
    '''
    sample_fname - a file with the first line 'ID' followed by one sample per line (7 digit number)
                  (no negatives or missings)
    '''
    all_samples = get_all_samples()

    with open(sample_fname) as samples_file:
        samples = np.array([line.strip() for line in samples_file][1:], dtype=int).reshape(-1, 1)

    samples_indicator = np.concatenate((samples, samples), axis=1)
    samples_merge = utils.merge_arrays(all_samples, samples_indicator)
    assert samples_merge.shape[1] == 2
    return samples_merge

def get_ordered_samples(sample_fname):
    arr = samples_array_with_indicator(sample_fname)
    return arr[~np.isnan(arr[:, 1]), 0]

def get_ordered_samples_ethnicity(ethnicity):
    return get_ordered_samples(f'{ukb}/sample_qc/runs/{ethnicity}/no_phenotype/combined_unrelated.sample')

def get_ordered_samples_phenotype(ethnicity, phenotype):
    return get_ordered_samples(f'{ukb}/sample_qc/runs/{ethnicity}/{phenotype}/combined_unrelated.sample')

def get_samples_idx(sample_fname):
    return ~np.isnan(samples_array_with_indicator(sample_fname)[:, 1])

def get_samples_idx_ethnicity(ethnicity):
    return get_samples_idx(f'{ukb}/sample_qc/runs/{ethnicity}/no_phenotype/combined_unrelated.sample')

def get_samples_idx_phenotype(ethnicity, phenotype):
    return get_samples_idx(f'{ukb}/sample_qc/runs/{ethnicity}/{phenotype}/combined_unrelated.sample')

