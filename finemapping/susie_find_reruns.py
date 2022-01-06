#!/usr/bin/env python3

import argparse
import collections
import json
import os.path

import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outfile')
# regions are tab sep
# cols chrom, start, end, any_strs
parser.add_argument('pheno_to_region_fname_json')
# dirs containing
# dirs for each run: {chr}_{start}_{end}
# each containing files:
# converged.txt containing the line TRUE or FALSE
# cs1.txt, ... , cs{n}.txt files  (n <= 50)
# no other files starting with 'cs'
# which have the members of the cs space sep on the first row
# and on the third row the first word is the min correlation between any members of the cs
parser.add_argument('pheno_to_susie_results_dir_json')
args = parser.parse_args()

pheno_to_region_fname = json.loads(args.pheno_to_region_fname_json)
pheno_to_susie_results_dir = json.loads(
    args.pheno_to_susie_results_dir_json
)

redo_regions = {}
for phenotype, region_fname in pheno_to_region_fname.items():
    regions = pl.read_csv(region_fname, sep='\t').to_dict(as_series=False)
    itr = zip(*(regions[col] for col in ('chrom', 'start', 'end', 'any_strs')))
    for chrom, start, end, any_strs in itr:
        if not any_strs:
            continue
        dir_ = (
            pheno_to_susie_results_dir[phenotype] +
            f'/{chrom}_{start}_{end}'
        )
        if not os.path.exists(dir_):
            redo_regions[(phenotype, chrom, start, end)] = 'missing dir'
            continue
        converged_fname = f'{dir_}/converged.txt'
        if not os.path.exists(converged_fname):
            redo_regions[(phenotype, chrom, start, end)] = 'missing converged'
            continue
        else:
            with open(converged_fname) as converged:
                did_converge = next(converged).strip()
                if did_converge == 'FALSE':
                    redo_regions[(phenotype, chrom, start, end)] = 'convergeance failed'
                else:
                    assert did_converge == 'TRUE'
        n_cses = len([fname for fname in os.listdir(dir_) if fname.startswith('cs')])
        if n_cses == 0:
            redo_regions[(phenotype, chrom, start, end)] = 'missing cses'
            continue
        if n_cses < 10:
            continue
        found_count = 0
        done = False
        redo = False
        for num in range(50):
            num+=1
            cs_fname = f'{dir_}/cs{num}.txt'
            if not os.path.exists(cs_fname):
                continue
            with open(cs_fname) as cs:
                n_vars = len(next(cs).split())
                next(cs)
                min_ld = float(next(cs).split()[0])
            if min_ld < .2 and n_vars > 50:
                done = True
                break
            found_count += 1
            if found_count == n_cses:
                done = True
                redo_regions[(phenotype, chrom, start, end)] = 'underexplored: no large nonsense CSes'
                break
        if not done:
            print(dir_)
            assert False

with open(args.outfile + '.debug', 'w') as outfile:
    for region, reason in redo_regions.items():
        outfile.write(f'{region} {reason}\n')
with open(args.outfile, 'w') as outfile:
    for (pheno, chrom, start, end) in redo_regions:
        outfile.write(
            pheno_to_susie_results_dir[pheno] +
            f'/{chrom}_{start}_{end}/alpha.tab\n'
        )
