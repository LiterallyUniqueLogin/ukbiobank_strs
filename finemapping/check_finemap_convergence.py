#!/usr/bin/env python3

import os

import polars as pl

import phenotypes

ukb = os.environ['UKB']

def generated_regions(phenotype):
    df = pl.scan_csv(
        f'{ukb}/signals/regions/{phenotype}.tab',
        sep='\t'
    ).filter('any_strs').collect()
    return [
        (chrom, start, end) for (chrom, start, end) in
        zip(df['chrom'], df['start'], df['end'])
        # remove a few for being too slow for finemapping
        if not ((
            phenotype == 'urate' and
            chrom == 4 and
            start == 8165642 and
            end == 11717761
        ) or (
            phenotype == 'total_bilirubin' and
            chrom == 12 and
            start == 19976272 and
            end == 22524428
        ))
    ]

for count, phenotype in enumerate(phenotypes.phenotypes_in_use):
    print(f"Working on phenotype #{count+1}: {phenotype}", flush=True)
    for (chrom, start, end) in generated_regions(phenotype):
        checked = False
        with open(f'{ukb}/finemapping/finemap_results/{phenotype}/{chrom}_{start}_{end}/finemap_output.log_sss') as logfile:
            for line in logfile:
                if 'converged after' in line:
                    n_iters = int(line.split()[-2])
                    assert n_iters < 100000 - 1
                    checked = True
                    break
        assert checked
