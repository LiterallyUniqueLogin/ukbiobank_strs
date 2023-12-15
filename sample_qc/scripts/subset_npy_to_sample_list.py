import argparse

import numpy as np

import python_array_utils
import sample_utils

parser = argparse.ArgumentParser()
parser.add_argument('npy')
parser.add_argument('sample_list')
args = parser.parse_args()

arr = np.load(args.npy)
samples = sample_utils.get_samples(args.sample_lsit)
arr = python_array_utils.merge_arrays(
    samples, arr
)

np.save('subsetted.npy', arr)
