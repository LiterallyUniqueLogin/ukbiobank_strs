import numpy as np

import python_array_utils as utils

def get_samples(samples_fname):
    with open(samples_fname) as samples_file:
        l1 = next(samples_file)
        samples = [line.split()[0] for line in samples_file]
    if l1.isdigit():
        samples = [l1] + samples
    samples = np.array(samples, dtype=int)
    samples = samples[samples != 0].reshape(-1, 1)
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

def get_samples_idx(all_samples_fname, samples_fname):
    return ~np.isnan(samples_array_with_indicator(all_samples_fname, samples_fname)[:, 1])

def n_samples(samples_fname):
    return get_samples(samples_fname).shape[0]

