#!/usr/bin/env python3

import argparse

import numpy as np

import python_array_utils as utils

parser = argparse.ArgumentParser()
parser.add_argument('outloc')
parser.add_argument('phenotype')
parser.add_argument('transformed_phenotype_file')
parser.add_argument('pheno_covar_names')
parser.add_argument('shared_covar_file')
parser.add_argument('shared_covar_names')
parser.add_argument('--conditional-genotypes')
parser.add_argument('--conditional-covar-names')
parser.add_argument('--binary', default=False, choices={'logistic', 'linear'})
args = parser.parse_args()

assert bool(args.conditional_genotypes) == bool(args.conditional_covar_names)

phenotype = args.phenotype

shared_covars = np.load(args.shared_covar_file)

subset_transformed_phenotype = np.load(args.transformed_phenotype_file)

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

with open(args.pheno_covar_names) as pheno_names_file:
    for line in pheno_names_file:
        line = line.strip()
        if line == '':
            continue
        col_names.append(line)

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

if not args.conditional_genotypes:
    np.savetxt(
        args.outloc,
        data,
        delimiter='\t',
        header='\t'.join(col_names),
        comments='',
        fmt='%.16g'
    )
else:
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

    np.savetxt(
        args.outloc,
        data,
        delimiter='\t',
        header='\t'.join(col_names),
        comments='',
        fmt='%.16g'
    )

