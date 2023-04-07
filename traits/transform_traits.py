#!/usr/bin/env python3

import argparse

import numpy as np
import scipy.stats

import python_array_utils as utils

def rank_phenotypes(data):
    """
    Ranks begin with 0
    """
    # "Ranking phenotype. Tie breaks for equal ranks are arbitrary and stable.\n")
    sort = np.argsort(data[:, 1])
    ranks = np.empty_like(sort)
    ranks[sort] = np.arange(data.shape[0])
    return ranks

def inverse_normalize_ranks(ranks):
    # "Inverse normalizing phenotype ranks to the standard normal distribution "
    # "via the transformation rank -> normal_quantile((rank + 0.5)/nsamples).\n")
    return scipy.stats.norm.ppf((ranks + 0.5)/ranks.shape[0])

def main():  # noqa: D103
    parser = argparse.ArgumentParser()
    parser.add_argument('outprefix')
    parser.add_argument('pheno_data')
    parser.add_argument('samples')
    parser.add_argument('--binary', default=False, action='store_true')
    args = parser.parse_args()

    # "Subsetting to samples with phenotype that passed sample_qc, "
    # f"as denoted by the file: {args.samples}\n")

    data = np.load(args.pheno_data)
    with open(args.samples) as sample_file:
        next(sample_file)
        samples = np.array([int(sample.strip()) for sample in sample_file])
    samples = samples.reshape(-1, 1)
    data = utils.merge_arrays(samples, data)


    # "Not standardizing covariates. Standardization should be done immediately before performing a regression. "
    # "This way I don't have to keep track of if covariates have standardized yet or not.\n"
    covariates = data[:, 2:]

    if not args.binary:
        ranks = rank_phenotypes(data)
        rin_ranks = inverse_normalize_ranks(ranks)
        transformed_data = np.concatenate(
            (samples, rin_ranks.reshape(-1, 1), covariates),
            axis=1
        )
    else:
        # "Binary outcome is left untransformed (0=case, 1=control)\n"
        phenotype = data[:, 1:2]
        transformed_data = np.concatenate(
            (samples, phenotype, covariates),
            axis=1
        )

    np.save(
        f'{args.outprefix}.npy',
        transformed_data
    )

if __name__ == "__main__":
    main()

