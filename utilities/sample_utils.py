import numpy as np

import python_array_utils as utils

def get_samples(samples_fname):
    with open(samples_fname) as samples_file:
        samples = np.array([line.split()[0] for line in samples_file][2:], dtype=int).reshape(-1, 1)
    return samples

def samples_array_with_indicator(all_samples_fname, samples_fname):
    all_samples = get_samples(all_samples_fname)
    samples = get_samples(samples_fname)

    samples_indicator = np.concatenate((samples, samples), axis=1)
    samples_merge = utils.merge_arrays(all_samples, samples_indicator)
    assert samples_merge.shape[1] == 2
    return samples_merge

def get_ordered_samples(all_samples_fname, samples_fname):
    arr = samples_array_with_indicator(all_samples_fname, samples_fname)
    return arr[~np.isnan(arr[:, 1]), 0]

'''
def get_ordered_samples_ethnicity(ethnicity):
    return get_ordered_samples(f'{ukb}/sample_qc/runs/{ethnicity}/no_phenotype/combined_unrelated.sample')

def get_ordered_samples_phenotype(ethnicity, phenotype):
    return get_ordered_samples(f'{ukb}/sample_qc/runs/{ethnicity}/{phenotype}/combined_unrelated.sample')
'''

def get_samples_idx(all_samples_fname, samples_fname):
    return ~np.isnan(samples_array_with_indicator(all_samples_fname, samples_fname)[:, 1])

def n_samples(samples_fname):
    return get_samples(samples_fname).shape[0]

'''
def get_samples_idx_ethnicity(ethnicity):
    return get_samples_idx(f'{ukb}/sample_qc/runs/{ethnicity}/no_phenotype/combined_unrelated.sample')

def get_samples_idx_phenotype(ethnicity, phenotype):
    return get_samples_idx(f'{ukb}/sample_qc/runs/{ethnicity}/{phenotype}/combined_unrelated.sample')
'''
