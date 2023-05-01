#!/usr/bin/env python3

import argparse

import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outdir')
parser.add_argument('phenotype')
parser.add_argument('my_str_gwas')
parser.add_argument('flank_start_to_start_and_end_pos')
parser.add_argument('hg19_pos_bed')
parser.add_argument('hg38_pos_bed')
parser.add_argument('repeat_units')
args = parser.parse_args()
phenotype = args.phenotype

flank_start_to_start_and_pos_table = pl.read_csv(
    args.flank_start_to_start_and_end_pos, sep='\t'
).groupby(
    ['chrom', 'pos', 'end_pos']
).agg([pl.col('snpstr_pos').first()]) # drop duplicate snpstr poses

hg19_pos_bed = pl.read_csv(
    args.hg19_pos_bed,
    sep='\t',
    has_header=False,
    new_columns=['chrom', 'pos', 'end_pos', 'ID']
).with_columns([
    pl.col('pos') + 1, # transform from bed to vcf coords
    pl.col('chrom').str.replace('chr', '').cast(int),
])
hg38_pos_bed = pl.read_csv(
    args.hg38_pos_bed,
    sep='\t',
    has_header=False,
    new_columns=['chrom', 'pos', 'end_pos', 'ID', 'drop']
).with_columns([
    pl.col('pos') + 1, # transform from bed to vcf coords
]).drop(['drop', 'chrom'])

hg19_pos_bed = hg19_pos_bed.join(
    hg38_pos_bed,
    how = 'left',
    on=['ID'],
    suffix='_hg38',
)

pos_table = hg19_pos_bed.join(
    flank_start_to_start_and_pos_table,
    on=['chrom', 'pos', 'end_pos']
).drop(['ID'])

repeat_units = pl.read_csv(args.repeat_units, sep='\t')

df = pl.read_csv(
    args.my_str_gwas, sep='\t'
).rename({'pos': 'snpstr_pos'}).join(
    pos_table,
    how='left',
    on=['chrom', 'snpstr_pos']
).join(
    repeat_units,
    how='left',
    on=['chrom', 'snpstr_pos']
).select([
    'chrom',
    pl.col('pos').alias('start_pos (hg19)'),
    pl.col('end_pos').alias('end_pos (hg19)'),
    pl.col('pos_hg38').alias('start_pos (hg38)'),
    pl.col('end_pos_hg38').alias('end_pos (hg38)'),
    'alleles',
    'locus_filtered',
    f'p_{phenotype}',
    f'coeff_{phenotype}', # todo is this and next on the wrong axis?
    f'se_{phenotype}',
    'R^2',
    pl.col('unit').alias('repeat_unit'), # use repeat unit as calculated in the paper instead of the trharmonizer way
    'period',
    'ref_len',
    # total / subset? - pref just subset
    # dosages sums, best guesses, best guess pairs? - pref all
    # het? entropy?
    # total dosage per summed gt
    # mean value per single dosgae
    # mean per paired dosage?
    pl.col('subset_HWEP').alias('hardy_weinberg_equilibrium_p_value'),
    'imputation_quality_per_allele_r2',
])

assert np.all(~df['start_pos (hg19)'].is_null().to_numpy())
assert np.all(~df['repeat_unit'].is_null().to_numpy())

df.write_csv(
    args.out, sep='\t'
)
