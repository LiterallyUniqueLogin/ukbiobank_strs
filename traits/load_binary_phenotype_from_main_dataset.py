#!/usr/bin/env python3

import argparse
import datetime
import os
import sys

import numpy as np
import pandas as pd

import python_array_utils as utils

ukb = os.environ['UKB']

def load_date_data_field(fname):
    date_array = utils.df_to_recarray(pd.read_csv(
        fname,
        header=0,
        delimiter='\t',
        names=['skip1', 'id', 'date', 'skip2'],
        index_col=False,
        parse_dates=['date'],
        infer_datetime_format=True
    )[['id', 'date']])
    # drop first and last rows which because there are some extra tabs in the file

    # calculate dates in days (pandas by default loads them in milliseconds)
    date_array = date_array.astype([('id', int), ('date', np.dtype('M8[D]'))])
    # same data, just as floats
    return np.stack((
        date_array['id'],
        date_array['date'].astype(float)
    ), axis=1)

parser = argparse.ArgumentParser()
parser.add_argument('phenotype_name')
parser.add_argument('phenotype_field_id')

args = parser.parse_args()

phenotype = args.phenotype_name

with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt', 'w') as unit_file:
    unit_file.write('binary_0=control_1=case\n')

with open(f'{ukb}/traits/phenotypes/{phenotype}_covar_names.txt', 'w') as covar_names:
    covar_names.write('age_at_diagnosis_days\n')

with open(f'{ukb}/traits/phenotypes/{phenotype}_README.txt', 'w') as readme:
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    data_fname = f'{ukb}/main_dataset/extracted_data/{phenotype}_{args.phenotype_field_id}.txt'
    readme.write(f"Run date: {today}\n")
    readme.write(
        f"Loading phenotype {phenotype} from txt file {data_fname} \n"
        "Text file is assumed to contain date first reported for all case samples, "
        "missing values for all control samples.\n"
        "Assumes coding https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=819 "
        "Dropping cases with dates 1900-01-01, 1901-01-01, 2037-07-07\n"
    )

    # cols: id, date first reported
    data = load_date_data_field(data_fname)

    drop_no_date = np.datetime64('1900-01-01', 'D').astype(float)
    drop_before_birth = np.datetime64('1901-01-01', 'D').astype(float)
    drop_future = np.datetime64('2037-03-03', 'D').astype(float)

    drop_no_date_idx = data[:, 1] == drop_no_date
    n_drop_no_date = np.sum(drop_no_date_idx)
    drop_before_birth_idx = data[:, 1] == drop_before_birth
    n_drop_before_birth = np.sum(drop_before_birth)
    drop_future_idx = data[:, 1] == drop_future
    n_drop_future = np.sum(drop_future)

    readme.write(
        f'Dropping {n_drop_no_date} cases without at a date (special value 1900-01-01)\n'
    )
    readme.write(
        f'Dropping {n_drop_before_birth} cases recorded at a date before the participant was '
        'born (special value 1901-01-01)\n'
    )
    readme.write(
        f'Dropping {n_drop_future} cases where the date is set in the future (special value '
        '2037-03-03)\n'
    )

    data[:, 1][drop_no_date_idx] = np.nan
    data[:, 1][drop_before_birth_idx] = np.nan
    data[:, 1][drop_future_idx] = np.nan

    is_case = ~np.isnan(data[:, 1])
    n_cases = np.sum(is_case)
    readme.write('Left with {n_cases} cases\n')

    # cols: id, is_case, date first reported
    data = np.concatenate((data, is_case.reshape(-1, 1)), axis=1)[:, [0,2,1]]

    year_of_birth = np.genfromtxt(
        f'{ukb}/main_dataset/extracted_data/year_of_birth_34.txt',
        delimiter='\t',
        skip_header=1
    )[:, 1:-1]
    month_of_birth = np.genfromtxt(
        f'{ukb}/main_dataset/extracted_data/month_of_birth_52.txt',
        delimiter='\t',
        skip_header=1
    )[:, 1:-1]

    # cols: id, is case, date first reported, year of birth, month of birth
    data = utils.merge_arrays(utils.merge_arrays(data, year_of_birth), month_of_birth)
    year_month_str_array = np.char.add(np.char.add(
        data[:, 3].astype(int).astype(str),
        '-'
    ),  np.char.rjust(data[:, 4].astype(int).astype(str), 2, fillchar='0')
    )
    print(year_month_str_array)
    birthdate_array = np.datetime64(year_month_str_array, 'D').astype(float)
    # cols: id, is case, date first reported, birth date
    data = np.concatenate((data, birthdate_array.reshape(-1, 1)), axis=1)[:, [0, 1, 2, 5]]

    is_birth_date = np.datetime64('1902-02-02', 'D').astype(float)
    first_year_date = np.datetime64('1903-03-03', 'D').astype(float)
    n_birth_date = np.sum(data[:, 2] ==  is_birth_date)
    n_first_year = np.sum(data[:, 2] == first_year_date)
    readme.write(
        f'Setting participant age for {n_birth_date} birth date records '
        '(special value 1902-02-02) to zero.\n'
    )
    readme.write(
        f'Setting participant age for {n_first_year} first year records '
        '(special value 1903-03-03) to zero.\n'
    )

    in_first_year = (data[:, 2] ==  is_birth_date) | (data[:, 2] == first_year_date)
    data[in_first_year, 2] = 0

    data[is_case & ~in_first_year, 2] = (
        data[is_case & ~in_first_year, 2] -
        data[is_case & ~in_first_year, 3]
    )

    readme.write(
        'Setting participant age for all other cases to diagnosis date (year, month, day) '
        'minus birth date (year, month, day=1), measured in number of days.\n'
    )
    # cols: id, is_case, age_at_diagnosis
    data = data[:, :2]

    np.save(f'{ukb}/traits/phenotypes/{phenotype}.npy', data)

