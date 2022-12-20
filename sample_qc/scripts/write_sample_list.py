#!/usr/bin/env python3

import argparse

import polars as pl

import python_array_utils

parser = argparse.ArgumentParser()
parser.add_argument('sc_fname')
parser.add_argument('outfname')
parser.add_argument('--value', default = 1, type=int)

args = parser.parse_args()

df = python_array_utils.load_extracted_data_as_pl(args.sc_fname)
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
).to_csv(
    args.outfname
)
