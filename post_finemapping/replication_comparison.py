#!/usr/bin/env python3

import os

import bokeh.io
import bokeh.layouts
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

associations_df = pl.concat([
    pl.scan_csv(
        f'{ukb}/post_finemapping/intermediate_results/finemapping_all_concordance_{phenotype}.tab',
        sep='\t',
        dtypes={
            **{f'{ethnicity}_p_val': float for ethnicity in other_ethnicities},
            **{f'{ethnicity}_coeff': float for ethnicity in other_ethnicities},
            **{f'{ethnicity}_se': float for ethnicity in other_ethnicities}
        }
    ).filter('is_STR') for phenotype in phenotypes.phenotypes_in_use
]).select([
    'phenotype',
    'chrom',
    'pos',
    'p_val',
    'coeff',
    'se',
    ((pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 0) & (pl.col('p_val') <= 5e-8)).alias('finemapped_susie'),
    ((pl.col('finemap_pip') >= 0.8) & (pl.col('p_val') <= 5e-8)).alias('finemapped_finemap'),
    *[f'{ethnicity}_p_val' for ethnicity in other_ethnicities],
    *[f'{ethnicity}_coeff' for ethnicity in other_ethnicities],
    *[f'{ethnicity}_se' for ethnicity in other_ethnicities],
]).collect()

df = associations_df.join(
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

shifted_figs = []
shared_x_axes = []
first_row=True
for name, condition in named_conditions[1:]:
    fig_row = []
    shifted_figs.append(fig_row)
    for eth_count, ethnicity in enumerate(other_ethnicities):
        base_results = df.filter(pl.col('p_val') > 1e-300).filter(~pl.col(f'{ethnicity}_coeff').is_nan()).select([
            (-pl.col('p_val').log10()).alias('x'),
            (-pl.col(f'{ethnicity}_p_val').log10()*(pl.col(f'{ethnicity}_coeff')/pl.col(f'{ethnicity}_coeff').abs())*(pl.col('coeff')/pl.col('coeff').abs())).alias('y')
        ])
        base_r2 = np.corrcoef(base_results['x'].to_numpy(), base_results['y'].to_numpy())[0,1]**2
        condition_results = df.filter(pl.col('p_val') > 1e-300).filter(~pl.col(f'{ethnicity}_coeff').is_nan()).filter(condition).select([
            (-pl.col('p_val').log10()).alias('x'),
            (-pl.col(f'{ethnicity}_p_val').log10()*(pl.col(f'{ethnicity}_coeff')/pl.col(f'{ethnicity}_coeff').abs())*(pl.col('coeff')/pl.col('coeff').abs())).alias('y')
        ])
        r2 = np.corrcoef(condition_results['x'].to_numpy(), condition_results['y'].to_numpy())[0,1]**2
        #print(f'{condition} r^2 for ethnicity {ethnicity}: {r2:.3}')

        shifted_r2s = []
        for i in range(-100, 101): #range(-10, 11):
            if i == 0:
                continue
            shifted_results = df.filter(pl.col('p_val') > 1e-300).filter(~pl.col(f'{ethnicity}_coeff').is_nan()).filter(condition.shift(i)).select([
                (-pl.col('p_val').log10()).alias('x'),
                (-pl.col(f'{ethnicity}_p_val').log10()*(pl.col(f'{ethnicity}_coeff')/pl.col(f'{ethnicity}_coeff').abs())*(pl.col('coeff')/pl.col('coeff').abs())).alias('y')
            ]).to_pandas()
            shifted_r2s.append(np.corrcoef(shifted_results['x'].to_numpy(), shifted_results['y'].to_numpy())[0,1]**2)
        fig_kws = dict(
            title=ethnicity,
            width=1000,
            height=800,
            x_axis_label = 'p-val r^2',
            toolbar_location = None
        )
        if not first_row:
            fig_kws['x_range'] = shared_x_axes[eth_count]
        fig = bokeh.plotting.figure(
            **fig_kws
        )
        if first_row:
            shared_x_axes.append(fig.x_range)
        base_quantile = np.sum(np.array(shifted_r2s) < base_r2)/len(shifted_r2s)
        quantile = np.sum(np.array(shifted_r2s) < r2)/len(shifted_r2s)
        fig.circle(shifted_r2s, np.zeros(len(shifted_r2s)))
        fig.line([base_r2, base_r2], [-1, 1], color='blue', legend_label=f'{base_quantile} quantile')
        fig.line([r2, r2], [-1, 1], color='red', legend_label=f'{quantile} quantile')
        fig.yaxis.visible=False
        fig_row.append(fig)
    first_row=False

        #print(f'Shifted {condition} r^2s for ethnicity {ethnicity}: ', [f'{r2:.3f}' for r2 in sorted(shifted_r2s)])
total_fig = bokeh.layouts.grid(shifted_figs)
bokeh.io.export_png(total_fig, filename=f'{ukb}/post_finemapping/results/continuous_replication_shifts.png')
bokeh.io.export_svg(total_fig, filename=f'{ukb}/post_finemapping/results/continuous_replication_shifts.svg')
exit()

for comparison_name, cond_wb, cond_rep in [
    ('p-val', -pl.col('p_val').log10(), lambda ethnicity: -pl.col(f'{ethnicity}_p_val').log10()),
    ('squared_z_score', (pl.col('coeff')/pl.col('se')).pow(2), lambda ethnicity: (pl.col(f'{ethnicity}_coeff')/pl.col(f'{ethnicity}_se')).pow(2)),
]:
    cols = []
    first = True
    first_row = True
    saved_y_ranges = []
    for (condition_name, condition) in named_conditions:
        row = []
        cols.append(row)
        for eth_count, ethnicity in enumerate(other_ethnicities):
            results = df.filter(pl.col('p_val') > 1e-300).filter(~pl.col(f'{ethnicity}_coeff').is_nan()).filter(condition).select([
                cond_wb.alias('x'),
                (cond_rep(ethnicity)*(pl.col(f'{ethnicity}_coeff')/pl.col(f'{ethnicity}_coeff').abs())*(pl.col('coeff')/pl.col('coeff').abs())).alias('y')
            ])
            r2 = np.corrcoef(results['x'].to_numpy(), results['y'].to_numpy())[0,1]**2
            fig_args = dict(
                title=f'{condition_name}, {ethnicity}, r^2 {r2:.3f}',
                y_axis_label=f'{ethnicity} {comparison_name}',
                x_axis_label=f'white-brit {comparison_name}',
                width=400,
                height=400,
            )
            if not first:
                fig_args['x_range'] = saved_x_range
            if not first_row:
                fig_args['y_range'] = saved_y_ranges[eth_count]
            fig = bokeh.plotting.figure(
                **fig_args
            )
            if first:
                first = False
                saved_x_range = fig.x_range
            if first_row:
                saved_y_ranges.append(fig.y_range)
            fig.circle(results['x'].to_numpy(), results['y'].to_numpy())
            fig.background_fill_color = None
            fig.border_fill_color = None
            fig.grid.grid_line_color = None
            fig.toolbar_location = None
            fig.title.text_font_size = '10px'
            fig.axis.axis_label_text_font_size = '20px'
            row.append(fig)
        first_row = False

    bokeh.io.export_png(bokeh.layouts.layout(cols), filename=f'{ukb}/post_finemapping/results/continuous_replication_{comparison_name}.png')
    bokeh.io.export_svg(bokeh.layouts.layout(cols), filename=f'{ukb}/post_finemapping/results/continuous_replication_{comparison_name}.svg')

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
