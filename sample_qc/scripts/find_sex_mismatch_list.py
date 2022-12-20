#!/usr/bin/env python3

import argparse

import polars as pl

import python_array_utils as utils

parser = argparse.ArgumentParser()
parser.add_argument('sc_genetic_sex')
parser.add_argument('sc_reported_sex')
parser.add_argument('outfname')

args = parser.parse_args()

rep = utils.load_extracted_data_as_pl(args.sc_reported_sex)
rep_col = rep.columns[1]
gen = utils.load_extracted_data_as_pl(args.sc_genetic_sex)
gen_col = gen.columns[1]

rep.filter(
    ~pl.col(rep_col).is_null()
).join(
    gen.filter(
        ~pl.col(gen_col).is_null()
    ),
    on=['eid']
).filter(
    pl.col(rep_col) != pl.col(gen_col)
).rename({
    'eid': 'ID'
}).select('ID').sort('ID').write_csv(args.outfname)
