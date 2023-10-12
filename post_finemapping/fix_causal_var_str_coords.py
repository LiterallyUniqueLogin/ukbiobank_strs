#!/usr/bin/env python3

import argparse

import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('chrom', type=int)
parser.add_argument('vars_and_betas')
parser.add_argument('flank_start_to_start_and_end_pos')
args = parser.parse_args()

fix_hipstr  = pl.read_csv(
    args.flank_start_to_start_and_end_pos,
    sep='\t'
)

with open(args.vars_and_betas) as vars_and_betas:
    vars_ = next(vars_and_betas).strip().split('\t')
    betas = next(vars_and_betas)

    done = False
    try:
        next(vars_and_betas)
    except StopIteration:
        done = True
    assert done

with open(args.vars_and_betas.split('/')[-1], 'w') as out:
    first = True
    for var in vars_:
        if not first:
            out.write('\t')
        first = False
        if 'STR' not in var:
            out.write(var)
        else:
            fix = fix_hipstr.filter((pl.col('snpstr_pos') == int(var.split('_')[1])) & (pl.col('chrom') == args.chrom))
            assert fix.shape[0] in {0, 1}
            if fix.shape[0] == 1:
                out.write(f'STR_{fix[0, "pos"]}')
            else:
                out.write(var)
    out.write('\n')
    out.write(betas)

