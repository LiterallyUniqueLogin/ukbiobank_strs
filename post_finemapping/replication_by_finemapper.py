#!/usr/bin/env python3

import argparse

import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('results_table')
args = parser.parse_args()

non_black_ethnicities = [ 'south_asian', 'chinese', 'irish', 'white_other']
ethnicities = ['black'] + non_black_ethnicities

df = pl.scan_csv(
    args.results_table,
    sep='\t'
).select([
    'association_p_value',
    'direction_of_association',
    'FINEMAP_pcausal',
    'SuSiE_CS_pcausal',
    'other_ethnic_association_ps',
    'other_ethnic_effect_directions'
]).filter(
    pl.col('association_p_value') > 1e-300
).with_columns([
    pl.col('other_ethnic_effect_directions').str.extract(
        '^([^,]+)', 1
    ).alias('black_effect_direction'),
    pl.col('other_ethnic_association_ps').str.extract(
        '^([^,]+)', 1
    ).cast(float).alias('black_p'),
    *[pl.col('other_ethnic_effect_directions').str.extract(
        '^([^ ]+ ){' + str(i+1) + '}([^,]+)', 2
    ).alias(f'{ethnicity}_effect_direction') for i, ethnicity in enumerate(non_black_ethnicities)],
    *[pl.col('other_ethnic_association_ps').str.extract(
        '^([^ ]+ ){' + str(i+1) + '}([^,]+)', 2
    ).cast(float).alias(f'{ethnicity}_p') for i, ethnicity in enumerate(non_black_ethnicities)],
]).with_columns([
     pl.when(
        pl.col('association_p_value') > 0
    ).then(
        -pl.col('association_p_value').log10()
    ).otherwise(
        300
    ).alias('-log_10_p'),
    *[pl.when(
        pl.col(f'{ethnicity}_p') > 0
    ).then(
        -pl.col(f'{ethnicity}_p').log10()
    ).otherwise(
        300
    ).alias(f'-{ethnicity}_log_10_p')
    for ethnicity in ethnicities]
]).with_columns([
    pl.when(
        pl.col('direction_of_association') == '-'
    ).then(
        -pl.col('-log_10_p')
    ).otherwise(
        pl.col('-log_10_p')
    ).alias('-log_10_p'),
    *[pl.when(
        pl.col(f'{ethnicity}_effect_direction') == 'NaN'
    ).then(
        0
    ).when(
        pl.col(f'{ethnicity}_effect_direction') == '-'
    ).then(
        -pl.col(f'-{ethnicity}_log_10_p')
    ).otherwise(
        pl.col(f'-{ethnicity}_log_10_p')
    ).alias(f'-{ethnicity}_log_10_p')
    for ethnicity in ethnicities]
]).collect()

for title, file_name_part, expr in [
    ('FINEMAP hits', 'finemap', pl.col('FINEMAP_pcausal') >= 0.8),
    ('SuSiE hits', 'susie', pl.col('SuSiE_CS_pcausal') >= 0.8),
    ('FINEMAP or SuSiE hits', 'finemap_or_susie', (pl.col('FINEMAP_pcausal') >= 0.8) | (pl.col('SuSiE_CS_pcausal') >= 0.8)),
    ('Hits in both FINEMAP and SuSiE', 'both_finemap_and_susie',  (pl.col('FINEMAP_pcausal') >= 0.8) & (pl.col('SuSiE_CS_pcausal') >= 0.8)),
    ('only FINEMAP hits', 'finemap_not_susie', (pl.col('FINEMAP_pcausal') >= 0.8) & (pl.col('SuSiE_CS_pcausal') < 0.8)),
    ('only SuSIE hits', 'susie_not_finemap', (pl.col('FINEMAP_pcausal') < 0.8) & (pl.col('SuSiE_CS_pcausal') >= 0.8)),
]:
    curr_df = df.filter(expr)
    figs = []
    for ethnicity in ethnicities:
        r2 = np.corrcoef(curr_df['-log_10_p'].to_numpy(), curr_df[f'-{ethnicity}_log_10_p'].to_numpy())[1,0]**2
        fig = bokeh.plotting.figure(
            width=400,
            height=400,
            tools='',
            x_axis_label='directional -log_10_white_brits_p',
            y_axis_label=f'directional -log_10_{ethnicity}_p',
            title=ethnicity if figs else title
        )
        if not figs:
            fig.add_layout(bokeh.models.Title(text=ethnicity), 'above')
        fig.add_layout(bokeh.models.Title(text=f'r^2: {r2:.2}'), 'below')
        if len(figs) == 4:
            fig.add_layout(bokeh.models.Title(
                text=f'All graphs exclude loci with white_brits_p_value <= 1e-330'
            ), 'below')
        fig.background_fill_color = None
        fig.border_fill_color = None
        fig.toolbar_location = None
        fig.circle(
            curr_df['-log_10_p'].to_numpy(), curr_df[f'-{ethnicity}_log_10_p'].to_numpy()
        )
        fig.line([-300, 300], [0,0], color='black')
        fig.line([0,0], [np.min(curr_df[f'-{ethnicity}_log_10_p'].to_numpy()), np.max(curr_df[f'-{ethnicity}_log_10_p'].to_numpy())], color='black')
        figs.append(fig)
    bokeh.io.export_png(bokeh.layouts.row(figs), filename=f'post_finemapping/results/replication_by_finemapper_{file_name_part}.png')

