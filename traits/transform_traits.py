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
    args = parser.parse_args()

    data = np.load(args.pheno_data)

    covariates = data[:, 2:]
    ranks = rank_phenotypes(data)
    rin_ranks = inverse_normalize_ranks(ranks)
    transformed_data = np.concatenate(
        (samples, rin_ranks.reshape(-1, 1), covariates),
        axis=1
    )

    np.save(
        f'{args.outprefix}.npy',
        transformed_data
    )

if __name__ == "__main__":
    main()

