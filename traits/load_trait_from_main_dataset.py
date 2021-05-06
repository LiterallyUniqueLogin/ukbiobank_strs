#!/usr/bin/env python3

import argparse
import datetime
import os
import subprocess as sp
import sys

import numpy as np

import python_array_utils as utils

# TODO what was reason some of the population was asked to come back and others
# were not? This affects interpretation of populations where phenotype value and
# age are not only from initail assessment

ukb = os.environ['UKB']

assessment_dict = {
    'init-assessment': 0,
    'repeat-assessment-1': 1,
    'imaging-visit': 2,
    'repeat-imaging-1': 3
}

parser = argparse.ArgumentParser()
parser.add_argument('phenotype_name')
parser.add_argument('unit')
parser.add_argument(
    'age',
    choices = set(assessment_dict.keys())
)
parser.add_argument(
    '--instance-id'
)

args = parser.parse_args()

assert args.unit != 'binary'

phenotype = args.phenotype_name

def enc_filename(data_request_id):
    return f'{ukb}/main_dataset/raw_data/fields{data_request_id}.ukb'

with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt', 'w') as unit_file:
    unit_file.write(f'{args.unit}\n')

with open(f'{ukb}/traits/phenotypes/{phenotype}_README.txt', 'w') as readme, \
        open(f'{ukb}/traits/phenotypes/{phenotype}_covar_names.txt', 'w') as covar_names:
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    age_csv_fname = f'{ukb}/main_dataset/extracted_data/assessment_age.csv'
    data_fname = f'{ukb}/main_dataset/extracted_data/{phenotype}.txt'
    readme.write(f"Run date: {today}\n")
    readme.write(
        "Loading phenotype {phenotype} from txt file "
        f" {data_fname} \n"
    )

    data = np.genfromtxt(
        data_fname,
        skip_header = 1,
        delimiter='\t'
    )
    # drop first and last rows which because of the way data is extracted
    # and then read by numpy are always nans
    data = data[:, 1:]
    data = data[:, :-1]

    covars = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')
    ages = covars[:, [0, -3, -2, -1]]
    if args.age is None:
        readme.write("Not using age as a covarite\n")
    elif args.age in assessment_dict:
        covar_names.write('age\n')
        readme.write(
            f"Choosing age for each participant corresponding to the visit "
            f"'{args.age}'. This is being loaded from the shared_covars file "
            f"{ukb}/traits/shared_covars/shared_covars.npy\n"
        )

        col = assessment_dict[args.age]
        if col >= 3:
            raise Exception("Doesn't support col == 3 right now")

        data = utils.merge_arrays(data, ages[:, [0, col + 1]])
        data = data[np.all(~np.isnan(data), axis=1), :]

    elif args.age == 'first-available':
        covar_names.write('age\n')
        readme.write(
            "Choosing age for each participant based on their age at the first "
            "visit for which this phenotype had a recorded value. If the participant "
            "had this phenotype measured at multiple visits, only the phenotype value and age "
            "at the first visit are being used. The age "
            "is being loaded from the shared_covars file "
            f"{ukb}/traits/shared_covars/shared_covars.npy\n"
        )
        with open(data_fname) as data_file:
            data_header = next(data_file)
        col_names = data_header.split()
        assert col_names[0] == 'eid'
        field_id = col_names[1].split('-')[0]
        assess_nums = []
        for idx, col_name in enumerate(col_names[1:]):
            assert len(col_name) == len(field_id) + 4
            assert col_name.startswith(field_id)
            assert col_name.endswith('.0')
            assess_num = int(col_name[-3])
            # temporarily skip assess_num == 3 because we need to refresh the data basket to get 
            # ages for that assessment
            if assess_num >= 3:
                assert idx != 0
                assert idx == len(col_names) - 2
                break
            assess_nums.append(assess_num)

        # going to rewrite the first column to be the correct phenotype value
        # going to rewrite the new last column to be the correct age
        # going to throw out the rest
        data = utils.merge_arrays(data, ages)
        data = utils.merge_arrays(data, np.full((data.shape[0], 1), np.nan))
        for col_num, assess_num in enumerate(assess_nums):
            nans = np.isnan(data[:, 1])
            data[nans, 1]  = data[nans, col_num + 1]
            data[nans, -1] = data[nans, col_num - 4]
        data = data[:, [0, 1, -1]]
    else:
        raise Exception("Age option not understood")

    np.save(data, f'{ukb}/traits/phenotypes/{phenotype}.npy')

