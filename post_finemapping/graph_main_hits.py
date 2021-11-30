#!/usr/bin/env python3

import argparse
import datetime
import json
import subprocess as sp

import bokeh.io
import bokeh.models
import bokeh.plotting
import polars

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile')
    #cols chrom, start_pos, phenotype
    parser.add_argument('hits_table')
    parser.add_argument('snpstr_correspondence_table') #cols chrom, pos, end_pos, snpstr_pos
    #parser.add_argument('per_pheno_summary_tables', help='json dict from pheno to summary table')
    parser.add_argument('per_pheno_results', help='json dict from pheno to results file')
    parser.add_argument('units', help='json dict from pheno to unit')

    args = parser.parse_args()

    # make sure same phenos everywhere
    assert set(args.units) == set(args.per_pheno_results)

    df = polars.read_csv(
        args.hits_table,
        sep='\t'
    )

    snpstr_correspondence = polars.read_csv(
        args.snpstr_corredpondence_table,
        sep='\t'
    )

    df = df.join(snpstr_correspondence, how='left', left_on=['chrom', 'start_pos'], right_on=['chrom', 'pos'])
    df['unit'] = 'NA'
    df['unit'] = df['unit'].cast(polars.Categorical)

    for phenotype in args.units:
        df[df['phenotype'] == phenotype, 'unit'] = args.units[phenotype]

    for idx in range(df.shape[0]):
        subprocess.

    units = json.loads(args.units)
    #summary_tables = json.loads(args.per_pheno_summary_tables)
    results = json.loads(args.per_pheno_results)



    today = datetime.datetime.now().strftime("%Y_%m_%d")
    manhattan_plot.add_layout(bokeh.models.Title(
        text=f"Plot creation date: {today}",
        align="right",
        text_font_size='18px'
    ))


if __name__ == '__main__':
    main()
