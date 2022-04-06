#!/usr/bin/env python3

import argparse
import os

import numpy.random

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('size', type=int)
parser.add_argument('rep')
args = parser.parse_args()

with open(f'{ukb}/sample_qc/runs/white_brits/no_phenotype/combined_unrelated.sample') as infile:
    lines = infile.readlines()

assert lines[0] == 'ID\n'

rng = numpy.random.default_rng()
with open(f'{ukb}/sample_qc/subpops/white_brits_subset_{args.size}_{args.rep}.txt', 'w') as out:
    out.write('ID\n')
    for idx in rng.choice(len(lines) - 1, size=args.size, replace=False):
        out.write(lines[idx+1])

