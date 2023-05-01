#!/usr/bin/env python3

import argparse

import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outprefix')
parser.add_argument('finemapping_regions')
parser.add_argument('str_loci')
args = parser.parse_args()

loci = pl.read_csv(
    args.str_loci, sep='\t', has_header=False
).unique().rename({'column_1': 'chrom', 'column_2': 'pos'})

for chrom in range(1, 23):
    dfs = [
        loci.filter(
            (pl.col('pos') >= start) & (pl.col('pos') <= end) & (pl.col('chrom') == chrom)
        )
        for (region_chrom, start, end) in zip(
            *pl.read_csv(
                args.finemapping_regions, sep='\t'
            ).select([
                'chrom', 'start', 'end']
            ).to_dict(False).values()
        ) if region_chrom == chrom
    ]
    if len(dfs) > 0:
        pl.concat(dfs).sort('pos').write_csv(f'{args.outfname}_chr{chrom}.tab', sep='\t')
    else:
        pl.DataFrame({'chrom': [], 'pos': []}).write_csv(f'{args.outfname}_chr{chrom}.tab', sep='\t')
