#!/usr/bin/env python3

import argparse

import numpy as np

import python_array_utils as utils

parser = argparse.ArgumentParser()
parser.add_argument('outloc')
parser.add_argument('phenotype')
parser.add_argument('transformed_phenotype_file')
parser.add_argument('pheno_covar_names')
parser.add_argument('shared_covar_file', nargs='?')
parser.add_argument('shared_covar_names', nargs='?')
parser.add_argument('--conditional-genotypes')
parser.add_argument('--conditional-covar-names')
parser.add_argument('--binary', default=False, choices={'logistic', 'linear'})
args = parser.parse_args()

assert bool(args.conditional_genotypes) == bool(args.conditional_covar_names)

phenotype = args.phenotype

subset_transformed_phenotype = np.load(args.transformed_phenotype_file)

if args.shared_covar_file:
    shared_covars = np.load(args.shared_covar_file)
    data = utils.merge_arrays(
        subset_transformed_phenotype,
        shared_covars
    )
else:
    data = subset_transformed_phenotype

data = data[np.all(~np.isnan(data), axis=1), :]
data = np.concatenate((data[:, 0:1], data), axis=1)

col_names = ['FID', 'IID']
if not args.binary:
    # phenotype has already been transformed, so no need to further modify
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

with open(args.pheno_covar_names) as pheno_names_file:
    for line in pheno_names_file:
        line = line.strip()
        if line == '':
            continue
        col_names.append(line)

if args.shared_covar_names:
    with open(args.shared_covar_names) as covar_names_file:
        for line in covar_names_file:
            line = line.strip()
            if line == '':
                continue
            col_names.append(line)

if not args.binary:
    suffix = ''
else:
    suffix = '_' + args.binary

if args.conditional_genotypes:
    genotypes = np.load(args.conditional_genotypes)
    data = utils.merge_arrays(data, genotypes)

    with open(args.conditional_covar_names) as conditional_names_file:
        line = next(conditional_names_file)
        first = True
        for word in line.split():
            if first:
                first = False
                continue
            col_names.append(word)

# standardize covariates for numerical stability
stds = data[:, 3:].std(axis=0)
stds[stds == 0] = 1 # simmply ignore covariates which are constant
data[:, 3:] = (data[:, 3:] - data[:, 3:].mean(axis=0))/stds

np.savetxt(
    args.outloc,
    data,
    delimiter='\t',
    header='\t'.join(col_names),
    comments='',
    fmt='%.16g'
)
