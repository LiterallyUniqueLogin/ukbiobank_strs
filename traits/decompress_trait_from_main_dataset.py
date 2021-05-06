#!/usr/bin/env python3

import argparse
import datetime
import os
import subprocess as sp
import sys

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('phenotype_name')
parser.add_argument('data_field_id')

args = parser.parse_args()

phenotype = args.phenotype_name

def enc_filename(data_request_id):
    return f'{ukb}/main_dataset/raw_data/fields{data_request_id}.ukb'

data_request_id = None
for data_request_id_attempt in {29170, 41414}:
    with open(enc_filename(data_request_id_attempt)) as fields_file:
        fields = ' '.join(fields_file.readlines()).replace('\n', ' ')
        if f' {args.data_field_id} ' in fields:
            data_request_id = data_request_id_attempt
            break
if data_request_id is None:
    raise ValueError(f"Couldn't find a data basked with data field id {args.data_field_id}")


out = sp.run(
    f'cd {ukb}/main_dataset/extracted_data || {{ echo "Can\'t cd to {ukb}/main_dataset/extracted_data" ; exit 1 ; }} ; '
    f'../../ukb_utilities/ukbconv ../raw_data/ukb{data_request_id}.enc_ukb txt '
    '-e../../ukb_utilities/encoding.ukb '
    f'-s{args.data_field_id} '
    f'-o{phenotype}',
    shell=True,
    capture_output=True
)

if out.returncode != 0:
    print(out.stdout.decode())
    print(out.stderr.decode(), file=sys.stderr)
    out.check_returncode()


    data = np.genfromtxt(
        f'{ukb}/main_dataset/extracted_data/{phenotype}.txt',
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

