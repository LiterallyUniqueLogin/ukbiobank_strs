import argparse

import numpy as np

import python_array_utils
import sample_utils

parser = argparse.ArgumentParser()
parser.add_argument('out')
parser.add_argument('npy')
parser.add_argument('sample_list')
args = parser.parse_args()

arr = np.load(args.npy)
assert not np.any(np.isnan(arr[:, 1])), 'error, nulls before subsetting'
samples = sample_utils.get_samples(args.sample_list)
arr = python_array_utils.merge_arrays(
    samples, arr
)
arr = arr[~np.isnan(arr[:, 1]), :]
assert not np.any(np.isnan(arr[:, 1])), 'error, nulls after subsetting'

np.save(args.out, arr)
