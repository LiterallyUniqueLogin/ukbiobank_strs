#!/usr/bin/env python3

import argparse

import polars as pl

import python_array_utils

parser = argparse.ArgumentParser()
parser.add_argument('outdir')
parser.add_argument(
    'white_brits_sample_fname',
    help='sep=" ", header, first col is ID'
)
parser.add_argument(
    'ethnicity_fname',
    help='header, sep="\\t", cols values: [eid, responses+], coding=https://biobank.ndph.ox.ac.uk/crystal/coding.cgi?id=1001'
)
args = parser.parse_args()

white_brits = pl.read_csv(args.white_brits_sample_fname, sep=' ', new_columns = ['id'])

ethnicities = python_array_utils.load_extracted_data_as_pl(args.ethnicity_fname)
#ethnicities = pl.read_csv(args.ethnicity_fname, sep='\t', dtype={'eid': int})
n_instances = ethnicities.shape[1] - 1
response_cols = [f'response{i}' for i in range(n_instances)]
ethnicities = ethnicities.rename(dict(zip(ethnicities.columns, ['id'] + response_cols)))

print(ethnicities.shape)

# remove people already listed as white brits
white_brits['marker'] = [True]*white_brits.shape[0]
ethnicities = ethnicities.join(white_brits, how='left', on=['id'])
ethnicities = ethnicities[ethnicities['marker'].is_null()]

print('post white brits')
print(ethnicities.shape)

for col in response_cols:
    ethnicities[ethnicities[col] == -3, col] = None
    ethnicities[ethnicities[col].is_in([4, 4001, 4002, 4003]), col] = 4 # African
    ethnicities[ethnicities[col].is_in([3001, 3002, 3003]), col] = 3000 # S Asian

# filter people who responded inconsistently
filters = []
for col1 in response_cols:
    for col2 in response_cols:
        if col1 == col2:
            continue

        filters.append(
            ethnicities[col1].is_null() |
            ethnicities[col2].is_null() |
            (ethnicities[col1] == ethnicities[col2])
        )

combo_filter = filters[0]
for filter_ in filters[1:]:
    combo_filter = combo_filter & filter_
ethnicities = ethnicities[combo_filter]

print('post inconsistent', ethnicities.shape)

# aggregate vals
ethnicities['val'] = ethnicities[response_cols[0]]
for col in response_cols[1:]:
    ethnicities = ethnicities.with_column(pl
        .when(ethnicities['val'].is_null())
        .then(ethnicities[col])
        .otherwise(ethnicities['val'])
        .alias('val')
    )

'''
for val in ethnicities[ethnicities['val'].is_null()]['id']:
    print(val)
'''
ethnicities = ethnicities[~ethnicities['val'].is_null()]

print('post null removal')
print(ethnicities.shape)

coding_names = {
    'irish': 1002,
    'white_other': 1003,
    'south_asian': 3000,
    'black': 4,
    'chinese': 5
}

ethnicities['ID_1'] = ethnicities['id']
ethnicities['ID_2'] = ethnicities['id']
for name, coding in coding_names.items():
    ethnicities[ethnicities['val'] == coding, ['ID_1', 'ID_2', 'missing', 'sex']].to_csv(
        f'{args.outdir}/{name}.sample', sep=' '
    )

