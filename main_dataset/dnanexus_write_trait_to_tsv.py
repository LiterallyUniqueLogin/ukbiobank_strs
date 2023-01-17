#!/usr/bin/env python3

import argparse
import os.path

# all the following from
# https://github.com/dnanexus/OpenBio/blob/master/UKB_notebooks/ukb-rap-pheno-basic.ipynb

import dxpy
import dxdata
#import pyspark

parser = argparse.ArgumentParser()
parser.add_argument('id', type=int)
parser.add_argument('dataset')
args = parser.parse_args()

#sc = pyspark.SparkContext()
#spark = pyspark.sql.SparkSession(sc)

dispensed_dataset = dxpy.find_one_data_object(
    typename='Dataset',
    name=os.path.basename(args.dataset),
    folder=os.path.dirname(args.dataset),
    #name='app*.dataset',
    #folder='/',
    name_mode='glob')
dispensed_dataset_id = dispensed_dataset['id']

dataset = dxdata.load_dataset(id=dispensed_dataset_id)
participant = dataset['participant']

# Returns all field objects for a given UKB showcase field id

def fields_for_id(field_id):
    from distutils.version import LooseVersion
    field_id = str(field_id)
    fields = participant.find_fields(name_regex=r'^p{}(_i\d+)?(_a\d+)?$'.format(field_id))
    return sorted(fields, key=lambda f: LooseVersion(f.name))

# Returns all field names for a given UKB showcase field id

def field_names_for_id(field_id):
    return [f.name for f in fields_for_id(field_id)]

df = participant.retrieve_fields(names=field_names_for_id(args.id), engine=dxdata.connect())

df.toPandas().to_csv('data.tsv', sep='\t', index=False)
