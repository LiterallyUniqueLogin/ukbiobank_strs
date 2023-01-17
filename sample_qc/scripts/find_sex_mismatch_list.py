#!/usr/bin/env python3

import argparse

import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('sc_genetic_sex')
parser.add_argument('sc_reported_sex')
parser.add_argument('outfname')

args = parser.parse_args()

rep = pl.read_csv(args.sc_reported_sex, sep='\t')
rep_col = rep.columns[1]
gen = pl.read_csv(args.sc_genetic_sex, sep='\t')
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
