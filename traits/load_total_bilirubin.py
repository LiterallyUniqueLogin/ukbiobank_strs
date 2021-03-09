#!/usr/bin/env python3

import datetime
import os

import numpy as np

import python_array_utils as utils

ukb = os.environ['UKB']

def load_total_bilirubin():
    """
    Load total bilirubin data for each sample

    Measured in umol/L

    Output array is 2D: id, total_bilirubin, age at measurement
    """
    with open(f'{ukb}/traits/phenotypes/total_bilirubin_unit.txt', 'w') as units:
        units.write('umol/L')

    with open(f'{ukb}/traits/phenotypes/total_bilirubin_README.txt', 'w') as readme, \
            open(f'{ukb}/traits/phenotypes/total_bilirubin_covar_names.txt', 'w') as covar_names:

        floc = f'{ukb}/main_dataset/extracted_data/bilirubin.csv'
        today = datetime.datetime.now().strftime("%Y_%m_%d")
        readme.write(f"Run date: {today}\n")
        readme.write(
            f"Adding total_bilirubin phenotype "
            f"and date of measurement, File: {floc}. "
            f"total_bilirubin is taken from the first assessment where it was "
            f"sampled, if any.\n"
        )
        readme.flush()
        # cols id", "tbil0", "tbil1"
        with open(floc) as bilirubin_csv:
            data = np.genfromtxt(
                (line.replace('"', '') for line in bilirubin_csv),
                skip_header=1,
                usecols=[0,3,4],
                delimiter=","
            )

        shared_covars = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')
        data = utils.merge_arrays(shared_covars[:, [0, -3, -2, -1]], data)

        # cols samples, bilirubin, age at meas
        out_array = np.concatenate(
            (data[:, 0:1], np.full((data.shape[0], 2), np.nan)), axis=1
        )

        use_tbil0 = ~np.isnan(data[:, -2])
        use_tbil1 = ~use_tbil0 & ~np.isnan(data[:, -1])
        out_array[use_tbil0, 1] = data[use_tbil0, -2]
        out_array[use_tbil1, 1] = data[use_tbil1, -1]

        # max bilirubin value right now is 144. I don't know enough to say that
        # this is too high, so no max filter
        # I've checked, none are negative

        readme.write("Adding age at bilirubin measurement as a covariate\n")
        readme.flush()
        out_array[use_tbil0, 2] = data[use_tbil0, 1]
        out_array[use_tbil1, 2] = data[use_tbil1, 2]
        covar_names.write('age\n')

        out_array = out_array[~np.any(np.isnan(out_array), axis=1), :]

        np.save(f'{ukb}/traits/phenotypes/total_bilirubin.npy', out_array)

if __name__ == "__main__":
    load_total_bilirubin()

