#!/usr/bin/env python3

#LD between rs7108857  and rs2155380

import argparse

import bgen_reader
import numpy as np

import sample_utils

parser = argparse.ArgumentParser()
parser.add_argument('bgen')
parser.add_argument('all_samples_file')
parser.add_argument('samples_file')
args = parser.parse_args()

samples = sample_utils.get_samples_idx(args.all_samples_file, args.samples_file)

bgen = bgen_reader.open_bgen(args.bgen, verbose=False, allow_complex=True)

dosages = []
for var in 119077003, 119080037:
    probs = np.squeeze(
        bgen.read(np.searchsorted(bgen.positions, var))
    )[samples, :]
    dosages.append(probs[:, 1] + 2*probs[:, 2])

print(np.corrcoef(dosages[0], dosages[1])[0, 1]**2)
