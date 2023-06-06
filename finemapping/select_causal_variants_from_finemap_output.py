#!/usr/bin/env python3

import argparse
import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outprefix')
parser.add_argument('finemap_output_log')
parser.add_argument('finemap_input_z')
parser.add_argument('creds', nargs='+')
args = parser.parse_args()

found_post_pr = 0
best_cred_num = None
best_prob = 0
for line in open(args.finemap_output_log):
    if found_post_pr == 0:
        if 'Post-Pr' in line:
            found_post_pr += 1
        continue
    elif found_post_pr == 1:
        found_post_pr += 1
        continue
    if line[0] == '-':
        break

    cred_num, _, prob = line.strip().split()
    prob = float(prob)
    if prob > best_prob:
        best_cred_num = cred_num
        best_prob = prob

assert best_cred_num is not None

one_cred = [cred for cred in args.creds if cred.split('.')[-1] == f'cred{best_cred_num}']
assert len(one_cred) == 1

with open(one_cred[0]) as cred:
    while True:
        line = next(cred)
        if line[0] == '#':
            continue
        if line[:5] == 'index':
            continue
        causal_vars = line.split()[1::2]
        break

input_z = pl.read_csv(args.finemap_input_z, sep=' ')

causal_betas = []
for causal_var in causal_vars:
    causal_betas.append(
        input_z.filter(pl.col('rsid') == causal_var)[0, 'beta']
    )

print('\t'.join(causal_vars))
print('\t'.join(str(beta) for beta in causal_betas))

