#!/usr/bin/env python3

import datetime
import os
import os.path

import numpy as np

import python_array_utils as utils

ukb = os.environ['UKB']

def load_covars():
    """
    Load the sex and population PC covariates. Loads ages at various
    assessments into a separate array

    Only columns that may be nan are the assessment ages

    Notes
    -----
    If interested in the source of this data, read the READMEs in the
    directories of the loaded files.
    """
    with open(f'{ukb}/traits/shared_covars/README.txt', 'w') as readme, \
            open(f'{ukb}/traits/shared_covars/covar_names.txt', 'w') as covar_names:
        today = datetime.datetime.now().strftime("%Y_%m_%d")
        readme.write(f"Run date: {today}\n")

        floc = f'{ukb}/microarray/ukb46122_cal_chr1_v2_s488282.fam'
        cols = (0, 4)
        readme.write(
            f"Loading participant ID index and sex covariate. File: {floc}, cols: {cols}\n"
            "Sex is encoded as 1 (male) or 2 (female)\n"
        )
        readme.flush()
        ids_and_sex = np.genfromtxt(
            floc,
            usecols=cols,
            delimiter=" "
        )
        readme.write(f"{ids_and_sex.shape[0]} total participants\n")
        covar_names.write('sex\n')

        floc = f'{ukb}/misc_data/EGA/ukb_sqc_v2.txt'
        cols = list(range(25, 65))

        [covar_names.write(f"pc{col}\n") for col in range(1, 41)]
        readme.write(
            f"Adding PC covariates 1-40. File: {floc}, cols: {cols}. Participants in same "
            f"order as previous file.\n"
        )
        readme.flush()
        pcs = np.genfromtxt(
            floc,
            delimiter=" ",
            usecols=cols
        )

        # these arrays are in the same row order, so just concatenate
        data = np.concatenate((ids_and_sex, pcs), axis=1)
        assert not np.any(np.isnan(data))

        # remove redacted samples which are denoted by negative id s
        readme.write("Removing redacted samples as indicated by negative id numbers.\n")
        data = data[data[:, 0] >= 0, :]
        # assert sex is either 1 (male) or 2 (female)
        assert np.all((data[:, 1] == 1) | (data[:, 1] == 2))

        readme.write('Standardizing covariates (subtracting mean, then dividing by standard deviation)\n')

        covars = data[:, 1:]
        covars = (covars - covars.mean(axis=0))/covars.std(axis=0)
        data[:, 1:] = covars

        np.save(f'{ukb}/traits/shared_covars/shared_covars.npy', data)

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
        age_file_name = f'{ukb}/main_dataset/extracted_data/assessment_ages_21003.txt'
        readme.write(
            f"Loading ages at assessments into assessment_ages.npy from file: {age_file_name}\n"
        )
        readme.flush()
        with open(age_file_name) as age_file:
            assessment_age = np.genfromtxt(
                age_file,
                skip_header=1,
                delimiter='\t'
            )[:, 1:-1]

        np.save(f'{ukb}/traits/shared_covars/assessment_ages.npy', assessment_age)


if __name__ == "__main__":
    load_covars()

