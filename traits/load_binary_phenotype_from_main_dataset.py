#!/usr/bin/env python3

import argparse
import datetime
import os

import numpy as np
import pandas as pd

import python_array_utils as utils

parser = argparse.ArgumentParser()
parser.add_argument('data_fname')
parser.add_argument('outprefix')
parser.add_argument('samples')
parser.add_argument('year_of_birth_fname')
parser.add_argument('month_of_birth_fname')
parser.add_argument('date_of_death_fname')
parser.add_argument('date_of_most_recent_first_occurrence_update')
parser.add_argument('--zero-one-neg-nan', action='store_true')

args = parser.parse_args()

outprefix = args.outprefix

def load_date_data_field(fname, extra_field = False):
    names=['skip1', 'id', 'date', 'skip2']
    if extra_field:
        names.append('skip3')
    date_array = utils.df_to_recarray(pd.read_csv(
        fname,
        header=0,
        delimiter='\t',
        names=names,
        index_col=False,
        parse_dates=['date'],
        infer_datetime_format=True
    )[['id', 'date']])
    # drop first and last rows which because there are some extra tabs in the file

    # calculate dates in days (pandas by default loads them in milliseconds)
    date_array = date_array.astype([('id', int), ('date', np.dtype('M8[D]'))])
    missing_times = np.isnan(date_array['date'])

    # same data, just as floats
    float_array = np.stack((
        date_array['id'],
        date_array['date'].astype(float)
    ), axis=1)
    # converting dates to floats doesn't preserve nans, so do it manually
    float_array[missing_times, 1] = np.nan
    return float_array

def load_0_1_neg_nan_field(fname):
    data = np.genfromtxt(
        fname,
        skip_header = 1,
        delimiter='\t'
    )[:, [1,2]]
    # don't include participants who didn't respond
    data = data[~np.isnan(data[:, 1]), :]
    # don't include participants who responded "don't know" or "won't tell"
    data = data[data[:, 1] >= 0, :]
    return data

with open(f'{outprefix}_unit.txt', 'w') as unit_file:
    unit_file.write('binary_0=control_1=case\n')

with open(f'{outprefix}_covar_names.txt', 'w') as covar_names:
    covar_names.write('current_age_or_age_at_death\n')

with open(f'{outprefix}_README.txt', 'w') as readme:
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    readme.write(f"Run date: {today}\n")
    readme.write(
        f"Loading from txt file {args.data_fname} \n"
        "Text file is assumed to contain date first reported for all case samples, "
        "missing values for all control samples.\n"
    )

    if not args.zero_one_neg_nan:
        # cols: id, date first reported
        data = load_date_data_field(args.data_fname)
        # cols: id, is_case
        data[:, 1] = ~np.isnan(data[:, 1])
    else:
        # cols: id, is_case
        data = load_0_1_neg_nan_field(args.data_fname)

    readme.write(f"Subsetting to samples at {args.samples}\n")
    # load samples that have passed qc
    samples = np.genfromtxt(
        args.samples,
        skip_header=1
    )

    # TODO test this
    # subset to those samples
    data = utils.merge_arrays(
        samples.reshape(-1, 1),
        data
    )

    year_of_birth = np.genfromtxt(
        args.year_of_birth_fname,
        delimiter='\t',
        skip_header=1
    )[:, 1:-1]
    month_of_birth = np.genfromtxt(
        args.month_of_birth_fname,
        delimiter='\t',
        skip_header=1
    )[:, 1:-1]

    # cols: id, is case, year of birth, month of birth
    data = utils.merge_arrays(utils.merge_arrays(data, year_of_birth), month_of_birth)
    missing_birth = np.isnan(data[:, 2]) | np.isnan(data[:, 3])
    assert np.sum(missing_birth) < 100
    data = data[~missing_birth, :]

    n_cases = np.sum(data[:, 1])
    readme.write(f'{n_cases} cases\n')

    year_month_str_array = np.char.add(np.char.add(
        data[:, 2].astype(int).astype(str),
        '-'
    ),  np.char.rjust(data[:, 3].astype(int).astype(str), 2, fillchar='0')
    )
    birthdate_array = year_month_str_array.astype('datetime64[D]').astype(float)
    # cols: id, is case, birth date in days since the epoch
    data = np.concatenate((data, birthdate_array.reshape(-1, 1)), axis=1)[:, [0, 1, 4]]

    date_of_death = load_date_data_field(
        args.date_of_death_fname,
        extra_field = True
    )
    # cols: id, is case, birth date in days since the epoch, date of death in days since the epoch
    data = utils.merge_arrays(data, date_of_death)

    # last day the first occurences were updated
    now = np.datetime64(args.date_of_most_recent_first_occurrence_update, 'D').astype(float)

    readme.write(
        f'Data last updated at {args.date_of_most_recent_first_occurrence_update}\n'
    )

    # date is not NaT and is less than or equal to the last time data was reported
    is_dead = data[:, 3] <= now

    # cols: id, is case, birth date in days since the epoch, age at death in days
    data[is_dead, 3] = data[is_dead, 3] - data[is_dead, 2]
    readme.write(
        'Setting participant age for participants who died before the last update to their '
        'death date (year, month, day) minus their birth date (year, month, day=1) measured '
        'in days.\n'
    )

    data[~is_dead, 2] = now - data[~is_dead, 2]
    data[is_dead, 2] = data[is_dead, 3]
    data = data[:, :3]
    # cols: id, is case, age in days (or age at death in days)

    readme.write(
        'Setting participant age for participants who were alive at the last update to '
        'that date minus their birth date (year, month, day=1), measured in days.\n'
    )

    np.save(f'{outprefix}.npy', data)

