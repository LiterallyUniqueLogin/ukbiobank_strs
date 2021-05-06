#!/usr/bin/env python3

import argparse
import datetime
import os
import subprocess as sp
import sys

import numpy as np

import python_array_utils as utils

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

phenotype = args.phenotype_name

def enc_filename(data_request_id):
    return f'{ukb}/main_dataset/raw_data/fields{data_request_id}.ukb'

with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt', 'w') as unit_file:
    unit_file.write(f'{args.unit}\n')

with open(f'{ukb}/traits/phenotypes/{phenotype}_README.txt', 'w') as readme, \
        open(f'{ukb}/traits/phenotypes/{phenotype}_covar_names.txt', 'w') as covar_names:
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    age_csv_fname = f'{ukb}/main_dataset/extracted_data/assessment_age.csv'
    data_file = f'{ukb}/main_dataset/extracted_data/{phenotype}.txt',
    readme.write(f"Run date: {today}\n")
    readme.write(
        "Loading phenotype {phenotype} from txt file "
        f" {data_file} \n"
    )

    data = np.genfromtxt(
        skip_header = 1,
        delimiter='\t'
    )

    covars = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')
    if args.age in assessment_dict:
        col = assessment_dict[args.age]
        if col >= 3:
            raise Exception("Doesn't support col == 3 right now")

        readme.write(
            f"Choosing age for each participant corresponding to the visit "
            f"'{args.age}'. This is being loaded from the shared_covars file "
            f"{ukb}/traits/shared_covars/shared_covars.npy\n"
        )
        data = utils.merge_arrays(data, covars[0, col - 3])
        data = data[np.all(~np.isnan(data), axis=1), :]

        covar_names.write('age\n')
    else:
        raise Exception("Age option not understood")

    np.save(data, f'{ukb}/traits/phenotypes/{phenotype}.npy')

