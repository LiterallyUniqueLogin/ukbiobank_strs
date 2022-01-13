#!/usr/bin/env python3

import argparse
import itertools
import numpy as np
import polars as pl

import h5py

parser = argparse.ArgumentParser()
parser.add_argument('hdf5_file')
parser.add_argument('var_names_file', help='space delimited, with header, column is "rsid" (even though names may not be rsids)')
parser.add_argument('chosen_var', type=int, help='0 indexed')
args = parser.parse_args()

chosen_var = args.chosen_var

var_names = pl.read_csv(
    args.var_names_file,
    sep=' '
).select('rsid').to_numpy().reshape(-1)

with h5py.File(args.hdf5_file) as gtsh5:
    gts = gtsh5['gts']
    assert gts.shape[1] > gts.shape[0]
    chosen_gts = gts[chosen_var, :]
    for var in itertools.chain(range(0, chosen_var), range(chosen_var+1, gts.shape[0])):
        if var % 100 == 99:
            print(f'Working on var {var+1} ... ', end='\r')
        prod = np.abs(np.corrcoef(chosen_gts, gts[var, :])[0, 1])
        if prod >= .8:
            print(f'{var_names[var]}: {prod}')
