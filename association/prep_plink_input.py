#!/usr/bin/env python3

import argparse
import os

import numpy as np

import python_array_utils as utils

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('phenotype')
parser.add_argument('--conditional')
parser.add_argument('--binary', default=False, choices={'logistic', 'linear'})
args = parser.parse_args()

phenotype = args.phenotype

shared_covars = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')

subset_transformed_phenotype = np.load(f'{ukb}/traits/subset_transformed_phenotypes/white_brits/{phenotype}.npy')

# shared covars here aren't necessarily going to have exactly mean 0 and std 1
# because they were standardized before subsetting, but that's okay
data = utils.merge_arrays(
    subset_transformed_phenotype,
    shared_covars
)
data = np.concatenate((data[:, 0:1], data), axis=1)

col_names = ['FID', 'IID']
if not args.binary:
    col_names.append(f'rin_{phenotype}')
else:
    col_names.append(phenotype)
    if args.binary == 'logistic':
        # plink expects and ouptuts a 1=control, 2=case encoding
        # instead of the 0=control, 1=case encoding we use elsewhere
        # (from: https://www.cog-genomics.org/plink/2.0/input
        # under the `--1` section)
        data[:, 2] += 1
    else:
        # standardize the binary trait so that
        # (a) plink treats it as a continuous trait and does a linear regression and
        # (b) plink doesn't have numeric instabilities
        data[:, 2] = (data[:, 2] - data[:,2].mean())/data[:, 2].std()

with open(f'{ukb}/traits/phenotypes/white_brits/{phenotype}_covar_names.txt') as pheno_names_file:
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

if not args.binary:
    suffix = ''
else:
    suffix = '_' + args.binary

if not args.conditional:
    np.savetxt(
        f'{ukb}/association/results/{phenotype}/plink_snp{suffix}/input/transformed_phenotype_and_covars.tab',
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
        f'{conditional_loc}_plink{suffix}.tab',
        data,
        delimiter='\t',
        header='\t'.join(col_names),
        comments='',
        fmt='%.16g'
    )

