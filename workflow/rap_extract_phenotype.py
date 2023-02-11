#!/usr/bin/env python3

import argparse
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument('ID', type=int)
parser.add_argument('--instances', type=int, help='number of instances, zero indexed')
parser.add_argument('--arrays', type=int, help='number of array elements, one indexed')

args = parser.parse_args()

field_names = [f'-ifield_names=p{args.ID}']
if args.instances:
    field_names = [f"{field_name}_i{instance}" for field_name in field_names for instance in range(args.instances + 1)]
if args.arrays:
    field_names = [f"{field_name}_a{array}" for field_name in field_names for array in range(1, args.arrays + 1)]

field_names_str = " ".join(field_names)

sp.run((
    'dx run table-exporter '
    '-idataset_or_cohort_or_dashboard=/app46122_20220823045256.dataset '
    '--folder imputed_strs_paper/main_dataset/extracted_data/ '
    f'-ioutput={args.ID} '
    '-ioutput_format=TSV '
    '-icoding_option=RAW '
    '-iheader_style=UKB-FORMAT '
    '-ientity=participant '
    '-ifield_names=eid '
    + field_names_str
),
    shell=True,
    check=True
)
