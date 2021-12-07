#!/usr/bin/env python3

import argparse
import ast
import json
import shutil

import numpy as np
import polars as pl

other_ethnicities=['black', 'southeast_asian', 'chinese', 'irish', 'white_other']

def reformat_dosage_dict_str(dict_str):
    d = ast.literal_eval(dict_str)
    d = {k: v for (k, v) in d.items() if not np.isnan(v) and v != 0}
    total_dosage = sum(d.values())
    #print(json.dumps({k: f'{v/total_dosage*100:.2f}%' for (k, v) in d.items()}).replace('"', "'"))
    return json.dumps({k: f'{v/total_dosage*100:.2f}%' for (k, v) in d.items()}).replace('"', "'")

def fix_cols(cols, phenotype):
    col_sightings = {}
    for col in cols:
        col = col.replace(phenotype, 'phenotype')
        if col not in col_sightings:
            col_sightings[col] = 1
            yield col
        else:
            yield col + '__' + str(col_sightings[col])
            col_sightings[col] += 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outtable')
    parser.add_argument('outreadme')
    parser.add_argument('pos_to_snpstr_pos')
    parser.add_argument('intable')
    parser.add_argument('inreadme')
    parser.add_argument('spot_test_fname_json_dict_fname')

    args = parser.parse_args()

    with open(args.spot_test_fname_json_dict_fname) as json_file:
        spot_test_fname_json_dict = next(json_file)

    shutil.copy(args.inreadme, args.outreadme)
    with open(args.outreadme, 'w+') as readme:
        readme.write(
            'other_ethnic_association_ps - association p-values for the other '
            'ethnicities in the order ' +
            ','.join(other_ethnicities) + '\n'
        )
        readme.write(
            'other_ethnic_effect_directions - direction of association (+/-) '
            'for the other ethnicities in the order ' +
            ','.join(other_ethnicities) +
            " (NaN if that ethnicity's p > 0.05)\n"
        )
        for ethnicity in other_ethnicities:
            readme.write(
                f'{ethnicity}_population_allele_frequencies - frequencies of each allele '
                "(by dosage) among the ethnicity's tested population\n"
            )

    hits = pl.scan_csv(args.intable, sep='\t').with_column(
        pl.col('white_brit_allele_frequencies').str.replace_all('"', "'")
    )
    cols = hits.columns

    hits = hits.join(
        pl.scan_csv(args.pos_to_snpstr_pos, sep='\t'),
        how='left',
        left_on=['chrom', 'start_pos'],
        right_on=['chrom', 'pos']
    )

    spot_tests_fnames = {
        tuple(key.split('__')): fname
        for key, fname in
        json.loads(spot_test_fname_json_dict).items()
    }

    spot_tests = {}
    for outer_ethnicity in other_ethnicities:
        spot_tests[outer_ethnicity] = pl.concat([
            (pl.scan_csv(
                    spot_test_fname,
                    sep='\t',
                    dtype={'alleles': str},
                    null_values=['nan'],
                    with_column_names=lambda cols: list(fix_cols(cols, phenotype))
                ).select([
                    pl.lit(phenotype).alias('phenotype'),
                    'chrom',
                    'pos',
                    pl.col('p_phenotype').cast(float).alias(f'{ethnicity}_p'),
                    pl.when(pl.col('p_phenotype') >= 0.05).then(np.nan).when(pl.col('coeff_phenotype') > 0).then(pl.lit('+')).otherwise(pl.lit('-')).alias(f'{ethnicity}_effect_direction'),
                    pl.col('subset_total_per_allele_dosages').apply(reformat_dosage_dict_str).alias(f'{ethnicity}_population_allele_frequencies')
                ]))
            for (phenotype, _, _, ethnicity), spot_test_fname
            in spot_tests_fnames.items()
            if ethnicity == outer_ethnicity
        ])

    for ethnicity in other_ethnicities:
        hits = hits.join(
            spot_tests[ethnicity],
            how='left',
            left_on=['phenotype', 'chrom', 'snpstr_pos'],
            right_on=['phenotype', 'chrom', 'pos']
        )

    hits = hits.with_columns([
        pl.sum([pl.col(f'{ethnicity}_p').cast(str) + pl.lit(', ') for ethnicity in other_ethnicities])
             .str.replace(', $', '').alias('other_ethnic_association_ps'),
        pl.sum([pl.col(f'{ethnicity}_effect_direction').cast(str) + pl.lit(', ') for ethnicity in other_ethnicities])
             .str.replace(', $', '').alias('other_ethnic_effect_directions')
    ])
    '''
    .with_columns([
        pl.col('subset_total_per_allele_dosages')
         .apply(reformat_dosage_dict_str)
         .alias('ethnic_allele_frequencies'),
        pl.when(pl.col('ethnic_coeff') > 0)
         .then(pl.lit('+'))
         .otherwise(pl.lit('-'))
         .alias('ethnic_effect_direction')
    ]).collect()

    index_cols=['phenotype', 'chrom', 'snpstr_pos']

    allele_freqs = (
        hits.map(lambda df:
            df.groupby(index_cols).pivot(pivot_column="ethnicity", values_column="ethnic_allele_frequencies").first()
        ).rename({col: f'{col}_population_allele_frequencies' for col in other_ethnicities})
    )

    effect_dirs = (
        hits.map(lambda df:
            df.groupby(index_cols).pivot(pivot_column="ethnicity", values_column="ethnic_effect_direction").first()
        ).select([
            index_cols,
            (pl.sum([pl.col(ethnicity).cast(str) + pl.lit(', ') for ethnicity in other_ethnicities])
             .str.replace(', $', '').alias('other_ethnic_effect_directions'))
        ])
    )

    assoc_ps = (
        hits.map(lambda df:
            df.groupby(index_cols).pivot(pivot_column="ethnicity", values_column="ethnic_p").first()
        ).select([
            index_cols,
            (pl.sum([pl.col(ethnicity).cast(str) + pl.lit(', ') for ethnicity in other_ethnicities])
             .str.replace(', $', '').alias('other_ethnic_association_ps'))
        ])
    )

    hits = (hits
        .groupby(index_cols)
        .join(allele_freqs, how='left', on=index_cols)
        .join(effect_dirs, how='left', on=index_cols)
        .join(assoc_ps, how='left', on=index_cols)
        .select([
            cols,
            'other_ethnic_association_ps',
            'other_ethnic_effect_directions',
            [ f'{col}_population_allele_frequencies' for col in cols]
        ])
        .collect()
    )
    '''

    hits = hits.select([
        *cols,
        'other_ethnic_association_ps',
        'other_ethnic_effect_directions',
        *[f'{ethnicity}_population_allele_frequencies' for ethnicity in other_ethnicities]
    ]).collect()
    assert hits.shape[0] == pl.read_csv(args.intable, sep='\t').shape[0]

    hits.to_csv(args.outtable, sep='\t',)

if __name__ == '__main__':
    main()
