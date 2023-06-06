#!/usr/bin/env python3

import argparse
import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outprefix')
parser.add_argument('colnames')
parser.add_argument('finemap_input_z')
parser.add_argument('CSes', nargs='+')
args = parser.parse_args()

all_vars = [line.strip() for line in open(args.colnames).readlines()]

causal_vars = []
causal_betas = []

for CS in args.CSes:
    lines = open(CS).readlines()
    assert len(lines) == 3
    if float(lines[2].split()[0]) < 0.8:
        continue
    causal_vars.append(all_vars[int(lines[0].split()[0]) - 1])

input_z = pl.read_csv(args.finemap_input_z, sep=' ')

for causal_var in causal_vars:
    causal_betas.append(
        input_z.filter(pl.col('rsid') == causal_var)[0, 'beta']
    )

print('\t'.join(causal_vars))
print('\t'.join(str(beta) for beta in causal_betas))

