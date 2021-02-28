#!/usr/bin/env python3

import os

import numpy as np

import python_array_utils as utils

ukb = os.environ['UKB']

def load_height():
    """
    Load height data for each sample
    Height is measured in cm.

    Writes out the following array to $UKB/traits/phenotypes/height.npy:
       2D float array of dimensions n_samples x 3. Columns are
       id, height, age at height measurement
    """
    with open(f'{ukb}/traits/phenotypes/height_unit.txt', 'w') as units:
        units.write('cm')

    with open(f'{ukb}/traits/phenotypes/height_README.txt', 'w') as readme:
        floc = f'{ukb}/main_dataset/extracted_data/height.txt'
        readme.write(
            f"Loading height phenotype and age at measurement. File: {floc}.\n"
            f"Only reporting first height measurement taken even if there "
            f"were measurements at multiple vists.\n"
        )
        readme.flush()

        # cols "id", "height", "height_sampling"
        height_data = np.genfromtxt(
            floc,
            skip_header=4,
            delimiter=" "
        )

        # taller than tallest person or shorter than shortest adult
        filter_extremes = np.logical_or(height_data[:, 1] > 274,
                                        height_data[:, 1] < 54)
        n_filtered = np.sum(filter_extremes)
        readme.write(f"Filtering {n_filtered} height values that are taller than "
                     "the world's tallest person or shorter than the shortest.\n")
        readme.flush()
        height_data[np.ix_(filter_extremes, [1,2])] = np.nan
        assert not np.any(np.isnan(height_data))

        shared_covars = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')

        height_data = utils.merge_arrays(shared_covars[:, [0, -3, -2, -1]], height_data)

        height_data = np.concatenate(
            (height_data, np.full((height_data.shape[0], 1), np.nan)),
            axis=1
        )
        assert height_data.shape[1] == 7

        start_idx = 1
        has_height = ~np.isnan(height_data[:, 4])
        height_data[has_height, 6] = height_data[
            has_height,
            start_idx + height_data[has_height, 5].astype(int)
        ]

        np.save(f'{ukb}/traits/phenotypes/height.npy', height_data[:, [0,4,6]])

if __name__ == "__main__":
    load_height()

