#!/usr/bin/env python3

import argparse
import datetime
import os
import subprocess as sp
import sys

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('field_name')
parser.add_argument('data_field_id')

args = parser.parse_args()

field_name = args.field_name
data_field_id = args.data_field_id

def enc_filename(data_request_id):
    return f'{ukb}/main_dataset/raw_data/fields{data_request_id}.ukb'

data_request_id = None
for data_request_id_attempt in {46781, 46782}:
    with open(enc_filename(data_request_id_attempt)) as fields_file:
        fields = ' '.join(fields_file.readlines()).replace('\n', ' ')
        if f' {data_field_id} ' in fields:
            data_request_id = data_request_id_attempt
            break
if data_request_id is None:
    raise ValueError(f"Couldn't find a data basket with data field id {data_field_id}")

with open(f'{ukb}/main_dataset/extracted_data/{field_name}_{data_field_id}_README.txt', 'w') as readme:
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    readme.write(
        f"Run date: {today}\n"
        f'Loading data field ID {data_field_id} from data request '
        f'{data_request_id} located at {ukb}/main_dataset/raw_data/ukb{data_request_id}.enc_ukb .'
        f"\nNaming this data field '{field_name}'. "
    )

out = sp.run(
    f'cd {ukb}/main_dataset/extracted_data || {{ echo "Can\'t cd to {ukb}/main_dataset/extracted_data" ; exit 1 ; }} ; '
    f'../../ukb_utilities/ukbconv ../raw_data/ukb{data_request_id}.enc_ukb txt '
    '-e../../ukb_utilities/encoding.ukb '
    f'-s{data_field_id} '
    f'-o{field_name}_{data_field_id}',
    shell=True,
    capture_output=True
)

if out.returncode != 0:
    print(out.stdout.decode())
    print(out.stderr.decode(), file=sys.stderr)
    out.check_returncode()


