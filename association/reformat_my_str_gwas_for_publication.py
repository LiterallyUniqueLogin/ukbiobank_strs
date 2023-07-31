#!/usr/bin/env python3

import argparse
import ast
import json

import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('out')
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

def to_frequencies(d_str):
    d = ast.literal_eval(d_str)
    total_dosage = sum(d.values())
    # strip braces
    return json.dumps({k: f'{v/total_dosage*100:.2f}%' for (k, v) in d.items()}).replace('"', "").replace(".0:", ":")[1:-1]

def to_number_of_common_alleles(d_str):
    d = ast.literal_eval(d_str)
    total_dosage = sum(d.values())
    return sum(v/total_dosage >= 0.01 for v in d.values())

def cleaned_dict_str(d_str):
    return d_str[1:-1].replace('"', '').replace('.0:', ':').replace('.0,', ',').replace('.0]', ']')

df = pl.read_csv(
    args.my_str_gwas,
    sep='\t',
    dtypes={'locus_filtered': str}
).filter(
    ((pl.col('chrom') != 17) | (pl.col('pos') != 80520458)) &
    ((pl.col('chrom') != 1) | (pl.col('pos') != 247747217)) &
    ((pl.col('chrom') != 1) | (pl.col('pos') != 247848392)) &
    ((pl.col('chrom') != 21) | (pl.col('pos') != 47741815)) &
    ((pl.col('chrom') != 8) | (pl.col('pos') != 145231731))
)

if 'subset_total_hardcall_alleles' in df.columns:
    df = df.rename({'subset_total_hardcall_alleles': 'subset_total_best_guess_alleles'})

if f'mean_{phenotype}_per_single_dosage' in df.columns:
    df = df.rename({
        f'mean_{phenotype}_per_single_dosage': f'mean_{phenotype}_per_summed_gt',
        '0.05_significance_CI': 'summed_0.05_significance_CI',
        '5e-8_significance_CI': 'summed_5e-8_significance_CI',
        f'mean_{phenotype}_per_paired_dosage': f'mean_{phenotype}_per_paired_gt',
        '0.05_significance_CI_duplicated_0': 'paired_0.05_significance_CI',
        '5e-8_significance_CI_duplicated_0': 'paired_5e-8_significance_CI',
    })

if all(df.filter(pl.col('chrom') == chrom)['pos'].is_in(pos_table.filter(pl.col('chrom') == chrom)['pos']).all() for chrom in range(1, 23)):
    df = df.join(
        pos_table,
        how='left',
        on=['chrom', 'pos']
    )
else:
    assert all(df.filter(pl.col('chrom') == chrom)['pos'].is_in(pos_table.filter(pl.col('chrom') == chrom)['snpstr_pos']).all() for chrom in range(1, 23))
    df = df.rename({'pos': 'snpstr_pos'}).join(
        pos_table,
        how='left',
        on=['chrom', 'snpstr_pos']
    )

df = df.join(
    repeat_units,
    how='left',
    on=['chrom', 'snpstr_pos']
).select([
    pl.col('chrom').alias('chromosome'),
    pl.col('pos').alias('base_pair_location'), # start_pos (hg19)
    pl.col('alleles').apply(lambda s: s.replace(".0", "")),
    pl.col(f'coeff_{phenotype}').alias('beta'),
    pl.col(f'se_{phenotype}').alias('standard_error'),
    pl.col('subset_total_per_allele_dosages').alias('allele_frequencies').apply(to_frequencies),
    pl.col(f'p_{phenotype}').alias('p_value'),
    'locus_filtered',
    #'imputation_quality_per_allele_r2', or maybe just correlation? alias info
    pl.col('ref_len').cast(str).alias('ref_allele').apply(lambda s: s.replace(".0", "")),
    pl.col('unit').alias('repeat_unit'), # use repeat unit as calculated in the paper instead of the trharmonizer way
    'period',
    pl.col('end_pos').alias('end_pos (hg19)'),
    pl.col('pos_hg38').alias('start_pos (hg38)'),
    pl.col('end_pos_hg38').alias('end_pos (hg38)'),
    pl.col('subset_total_best_guess_alleles').alias('n').apply(lambda d: sum(ast.literal_eval(d).values())//2),
    pl.col('subset_total_per_allele_dosages').alias('number_of_common_alleles').apply(to_number_of_common_alleles),
    pl.col(f'mean_{phenotype}_per_summed_gt').apply(cleaned_dict_str),
    pl.col('summed_0.05_significance_CI').apply(cleaned_dict_str),
    pl.col('summed_5e-8_significance_CI').apply(cleaned_dict_str),
    pl.col(f'mean_{phenotype}_per_paired_gt').apply(cleaned_dict_str),
    pl.col('paired_0.05_significance_CI').apply(cleaned_dict_str),
    pl.col('paired_5e-8_significance_CI').apply(cleaned_dict_str)
    # not including: model R^2, per allele r2, all sample pop values (just tested samples), single or paired best guesses (just single dosages)
])

assert np.all(~df['base_pair_location'].is_null().to_numpy())
assert np.all(~df['repeat_unit'].is_null().to_numpy())

unfiltered = df.filter(
    pl.col('locus_filtered') == 'False'
).drop('locus_filtered')

filtered = df.filter(
    pl.col('locus_filtered') != 'False'
)

unfiltered.write_csv(
    f'{args.out}.tab', sep='\t'
)
filtered.write_csv(
    f'{args.out}.filtered.tab', sep='\t'
)
