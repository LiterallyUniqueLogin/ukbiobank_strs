#!/usr/bin/env python3

import os

import bokeh.io
import bokeh.models
import bokeh.plotting
import numpy as np
import polars as pl
import scipy.stats

import phenotypes

ukb = os.environ['UKB']

other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']

causal_STR_candidates = pl.read_csv(
    f'{ukb}/post_finemapping/intermediate_results/concordant_causal_STR_candidates.tab',
    sep='\t'
).select([
    'phenotype',
    'chrom',
    'pos',
    pl.lit(True).alias('causal_STR_candidate_indicator')
])

p_vals_df = pl.concat([
    pl.scan_csv(
        f'{ukb}/post_finemapping/intermediate_results/finemapping_all_concordance_{phenotype}.tab',
        sep='\t',
        dtypes={f'{ethnicity}_p_val': float for ethnicity in other_ethnicities}
    ).filter('is_STR') for phenotype in phenotypes.phenotypes_in_use
]).select([
    'phenotype',
    'chrom',
    'pos',
    'p_val',
    ((pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 0) & (pl.col('p_val') <= 5e-8)).alias('finemapped_susie'),
    ((pl.col('finemap_pip') >= 0.8) & (pl.col('p_val') <= 5e-8)).alias('finemapped_finemap'),
    *[f'{ethnicity}_p_val' for ethnicity in other_ethnicities],
]).collect()

df = p_vals_df.join(
    causal_STR_candidates,
    how='left',
    on=['phenotype', 'chrom', 'pos']
).with_column((~pl.col('causal_STR_candidate_indicator').is_null()).alias('is_causal_STR_candidate'))


named_conditions = [
    ('gwsig STRs', pl.col('p_val') <= 5e-8),
    ('SuSiE fine-mapped STRs', pl.col('finemapped_susie')),
    ('FINEMAP fine-mapped STRs', pl.col('finemapped_finemap')),
    ('causal STR candidates', pl.col('is_causal_STR_candidate'))
]

cols = []
for (name, condition) in named_conditions:
    row = []
    cols.append(row)
    for ethnicity in other_ethnicities:
        results = df.filter(pl.col('p_val') > 1e-300).filter(condition).select([(-pl.col('p_val').log10()).alias('x'), (-pl.col(f'{ethnicity}_p_val').log10()).alias('y')])
        r2 = np.corrcoef(results['x'].to_numpy(), results['y'].to_numpy())[0,1]**2
        fig = bokeh.plotting.figure(
            title=f'{name}, {ethnicity}, r^2 {r2:.3f}',
            y_axis_label=f'{ethnicity} p-value',
            x_axis_label='white-brit p-value',
            width=400,
            height=400,
        )
        fig.circle(results['x'].to_numpy(), results['y'].to_numpy())
        fig.background_fill_color = None
        fig.border_fill_color = None
        fig.grid.grid_line_color = None
        fig.toolbar_location = None
        fig.title.text_font_size = '10px'
        fig.axis.axis_label_text_font_size = '20px'
        row.append(fig)
    bokeh.io.export_png(bokeh.layouts.layout(cols), filename=f'{ukb}/post_finemapping/results/continuous_replication.png')
    bokeh.io.export_svg(bokeh.layouts.layout(cols), filename=f'{ukb}/post_finemapping/results/continuous_replication.svg')

'''
for ethnicity in other_ethnicities:
    for name, condition in named_conditions
    n_neither = df.filter(pl.col('p_val') <= 5e-8).filter((pl.col(f'{ethnicity}_p_val') > 0.05) & ~pl.col('is_causal_STR_candidate')).shape[0]
    n_replicate_not_candidates = df.filter(pl.col('p_val') <= 5e-8).filter((pl.col(f'{ethnicity}_p_val') <= 0.05) & ~pl.col('is_causal_STR_candidate')).shape[0]
    n_candidates_not_replicate = df.filter(pl.col('p_val') <= 5e-8).filter((pl.col(f'{ethnicity}_p_val') > 0.05) & pl.col('is_causal_STR_candidate')).shape[0]
    n_both = df.filter(pl.col('p_val') <= 5e-8).filter((pl.col(f'{ethnicity}_p_val') <= 0.05) & pl.col('is_causal_STR_candidate')).shape[0]

    contingency_table = [[n_neither, n_replicate_not_candidates], [n_candidates_not_replicate, n_both]]
    if np.any(np.array(contingency_table) < 20):
        p_val = scipy.stats.fisher_exact(contingency_table)[1]
    else:
        p_val = scipy.stats.chi2_contingency(contingency_table)[1]
    print(ethnicity, p_val)

xs = []
ys = []
for ethnicity in other_ethnicities:
    for (name, condition) in named_conditions:
        xs.append((ethnicity, name))
        ys.append(
            df.filter(
                condition
            ).select(
                (pl.col(f'{ethnicity}_p_val') <= 0.05).cast(int).mean().alias('out')
            )['out'].to_numpy()[0]
        )

fig = bokeh.plotting.figure(
    title='Replication rate by category and ethnicity',
    x_range=bokeh.models.FactorRange(*xs),
    y_axis_label = 'Replication rate',
    width=1200,
    height=800,
)
fig.vbar(x=xs, top=ys, width=0.9)
fig.y_range.start=0
fig.x_range.range_padding=0.1
fig.xaxis.major_label_orientation = 1
fig.background_fill_color = None
fig.border_fill_color = None
fig.grid.grid_line_color = None
fig.toolbar_location = None
fig.title.text_font_size = '30px'
fig.axis.axis_label_text_font_size = '26px'
fig.axis.major_label_text_font_size = '20px'
bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/binary_replication.png')
bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/binary_replication.svg')
'''
