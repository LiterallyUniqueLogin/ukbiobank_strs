#!/usr/bin/env python3

import argparse

import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('sc_fname')
parser.add_argument('outfname')
parser.add_argument('--value', default = 1, type=int)

args = parser.parse_args()

df = pl.read_csv(args.sc_fname, sep='\t')
print(df.columns[1])
print(df.filter(
    ~pl.col(df.columns[1]).is_null()
))
print(df.with_column(
    pl.col(df.columns[1]).cast(int)
).filter(
    ~pl.col(df.columns[1]).is_null()
))
df.with_column(
    pl.col(df.columns[1]).cast(int)
).filter(
    ~pl.col(df.columns[1]).is_null() &
    (pl.col(df.columns[1]) == args.value)
).select(
    'eid'
).rename(
    {'eid': 'ID'}
).sort(
    'ID'
).write_csv(
    args.outfname
)
