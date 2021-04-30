#!/usr/bin/env python3

import argparse
import os

import numpy as np

import python_array_utils as utils

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('phenotype')
parser.add_argument('--conditional')
args = parser.parse_args()

phenotype = args.phenotype

shared_covars = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')[:, :-3]

subset_rin_phenotype = np.load(f'{ukb}/traits/subset_rin_phenotypes/{phenotype}.npy')

data = utils.merge_arrays(
    subset_rin_phenotype,
    shared_covars
)
data = np.concatenate((data[:, 0:1], data), axis=1)

col_names = ['FID', 'IID', f'rin_{phenotype}']
with open(f'{ukb}/traits/phenotypes/{phenotype}_covar_names.txt') as pheno_names_file:
    for line in pheno_names_file:
        line = line.strip()
        if line == '':
            continue
        col_names.append(line)

with open(f'{ukb}/traits/shared_covars/covar_names.txt') as covar_names_file:
    for line in covar_names_file:
        line = line.strip()
        if line == '':
            continue
        col_names.append(line)

if not args.conditional:
    np.savetxt(
        f'{ukb}/association/results/{phenotype}/plink_snp/input/rin_phenotype_and_covars.tab',
        data,
        delimiter='\t',
        header='\t'.join(col_names),
        comments='',
        fmt='%.16g'
    )
else:
    conditional_loc = f'{ukb}/association/results/{args.phenotype}/conditional_inputs/{args.conditional}'
    genotypes = np.load(f'{conditional_loc}.npy')
    data = utils.merge_arrays(data, genotypes)

    with open(f'{conditional_loc}_varnames.txt') as conditional_names_file:
        line = next(conditional_names_file)
        first = True
        for word in line.split():
            if first:
                first = False
                continue
            col_names.append(word)

    np.savetxt(
        f'{conditional_loc}_plink.txt',
        data,
        delimiter='\t',
        header='\t'.join(col_names),
        comments='',
        fmt='%.16g'
    )
