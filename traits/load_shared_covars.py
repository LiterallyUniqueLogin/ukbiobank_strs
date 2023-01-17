#!/usr/bin/env python3

import argparse

import numpy as np

import python_array_utils

def load_covars():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('fam_file')
    parser.add_argument('pcs_fname')
    parser.add_argument('assessment_ages_file')
    args = parser.parse_args()

    """
    Load the sex and population PC covariates. Loads ages at various
    assessments into a separate array

    Only columns that may be nan are the assessment ages
    """
    with open(f'{args.outdir}/covar_names.txt', 'w') as covar_names:
        cols = (0, 4)
        # Sex is encoded as 1 (male) or 2 (female)
        ids_and_sex = np.genfromtxt(
            args.fam_file,
            usecols=cols,
            delimiter=" "
        )
        covar_names.write('sex\n')
        # remove redacted samples which are denoted by negative id s
        ids_and_sex = ids_and_sex[ids_and_sex[:, 0] > 0, :]
        assert not np.any(np.isnan(ids_and_sex))

        pcs = np.genfromtxt(
            args.pcs_fname,
            delimiter="\t",
            skip_header = 1,
        )
        [covar_names.write(f"pc{col}\n") for col in range(1, 41)]

        data = python_array_utils.merge_arrays(
            ids_and_sex, pcs
        )
        print(data.shape)
        print('IDs without PCs (that are being skipped)')
        print(list(data[np.any(np.isnan(data), axis=1), 0]))
        data = data[~np.any(np.isnan(data), axis=1), :]
        print(data.shape)

        # assert sex is either 1 (male) or 2 (female)
        assert np.all((data[:, 1] == 1) | (data[:, 1] == 2))

        # TODO stop this
        # Standardizing covariates (subtracting mean, then dividing by standard deviation)
        covars = data[:, 1:]
        covars = (covars - covars.mean(axis=0))/covars.std(axis=0)
        data[:, 1:] = covars

        np.save(f'{args.outdir}/shared_covars.npy', data)

        # Age will be included as a covariate for all dependent variables
        # to prevent confounding.
        # If not, then if age was correlated with the dependent variable,
        # any loci which caused people to live longer or shorter would be
        # spuriously associated with the dependent variable.
        # However, each dependent variable measured during an assessment may
        # have been measured at one of multiple assessments.
        # As such, each measured dependent variable needs to specify which assessments
        # the age should be drawn from for each participant.
        # because of that, ages are saved as a separate file
        age_file_name = args.assessment_ages_file
        with open(age_file_name) as age_file:
            assessment_age = np.genfromtxt(
                age_file,
                skip_header=1,
                delimiter='\t'
            )

        np.save(f'{args.outdir}/assessment_ages.npy', assessment_age)

if __name__ == "__main__":
    load_covars()

