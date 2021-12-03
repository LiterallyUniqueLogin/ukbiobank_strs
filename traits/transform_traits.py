#!/usr/bin/env python3

import argparse
import datetime
import os
import os.path

import numpy as np
import scipy.stats

import python_array_utils as utils

ukb = os.environ['UKB']

def rank_phenotypes(readme, data):
    """
    Ranks begin with 0
    """
    readme.write("Ranking phenotype. Tie breaks for equal ranks are arbitrary and stable.\n")
    readme.flush()
    sort = np.argsort(data[:, 1])
    ranks = np.empty_like(sort)
    ranks[sort] = np.arange(data.shape[0])
    return ranks

def inverse_normalize_ranks(readme, ranks):
    readme.write("Inverse normalizing phenotype ranks to the standard normal distribution "
                 "via the transformation rank -> normal_quantile((rank + 0.5)/nsamples).\n")
    readme.flush()
    return scipy.stats.norm.ppf((ranks + 0.5)/ranks.shape[0])

def main():  # noqa: D103
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('ethnicity')
    parser.add_argument('--binary', default=False, action='store_true')
    args = parser.parse_args()

    with open(f'{ukb}/traits/subset_transformed_phenotypes/{args.ethnicity}/{args.phenotype}_README.txt', 'w') as readme:
        today = datetime.datetime.now().strftime("%Y_%m_%d")
        readme.write(f"Run date: {today}\n")
        sample_fname = f'{ukb}/sample_qc/runs/{args.ethnicity}/{args.phenotype}/combined_unrelated.sample'
        readme.write("Subsetting to samples with phenotype that passed sample_qc, "
                     f"as denoted by the file: {sample_fname}\n")
        readme.flush()

        data = np.load(f'{ukb}/traits/phenotypes/{args.ethnicity}/{args.phenotype}.npy')
        with open(sample_fname) as sample_file:
            next(sample_file)
            samples = np.array([int(sample.strip()) for sample in sample_file])
        samples = samples.reshape(-1, 1)
        data = utils.merge_arrays(samples, data)


        readme.write(
            "Standardizing covariates (subtracting mean, then dividing by standard deviation)\n"
        )
        readme.flush()
        covariates = data[:, 2:]
        standardized_covariates = \
                (covariates - covariates.mean(axis=0))/covariates.std(axis=0)

        if not args.binary:
            ranks = rank_phenotypes(readme, data)
            rin_ranks = inverse_normalize_ranks(readme, ranks)
            transformed_data = np.concatenate(
                (samples, rin_ranks.reshape(-1, 1), standardized_covariates),
                axis=1
            )
        else:
            readme.write(
                "Binary outcome is left untransformed (0=case, 1=control)\n"
            )
            phenotype = data[:, 1:2]
            transformed_data = np.concatenate(
                (samples, phenotype, standardized_covariates),
                axis=1
            )

        np.save(
            f'{ukb}/traits/subset_transformed_phenotypes/{args.ethnicity}/{args.phenotype}.npy',
            transformed_data
        )

if __name__ == "__main__":
    main()

