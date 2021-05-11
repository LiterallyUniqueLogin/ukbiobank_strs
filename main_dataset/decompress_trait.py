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

with open(f'{ukb}/main_dataset/extracted_data/{phenotype}_README.txt', 'w') as readme:
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    readme.write(
        f"Run date: {today}\n"
        f'Loading data field ID {args.data_field_id} from data request '
        f'{data_request_id} located at {ukb}/main_dataset/raw_data/ukb{data_request_id}.enc_ukb .'
        f"\nNaming this data field '{phenotype}'. "
    )

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


