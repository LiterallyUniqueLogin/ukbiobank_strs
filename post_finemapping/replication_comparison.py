#!/usr/bin/env python3

import os

import bokeh.embed
import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import bokeh.resources
import bokeh.transform
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import scipy.stats
import statsmodels.formula.api
from statsmodels.stats.proportion import proportions_ztest

import phenotypes

ukb = os.environ['UKB']

other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']
other_ethnicities_w_counts = ['black\nn=7,562', 'south_asian\nn=7,397', 'chinese\nn=1,525', 'irish\nn=11,978', 'white_other\nn=15,838']

# from https://stackoverflow.com/a/45846841
def human_format(num):
    num = float('{:.2g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'k', 'm', 'b', 't'][magnitude])

print('loading data', flush=True)

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
    'region',
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
    ('genome-wide significant STRs', pl.col('p_val') <= 5e-8),
    ('FINEMAP STRs with PIP >= 0.8', pl.col('finemapped_finemap')),
    ('SuSiE STRs with max alpha >= 0.8', pl.col('finemapped_susie')),
    ('confidently fine-mapped STRs', pl.col('is_causal_STR_candidate'))
]

print("performing logistic regressions", flush=True)
for ethnicity in other_ethnicities:
    for group_comp in ((0, 1), (0, 2), (0, 3), (1, 3), (2, 3)):
        category_df = df.filter(
            (pl.col('p_val') <= 1e-10) &
            named_conditions[group_comp[0]][1] &
            ~pl.col(f'{ethnicity}_coeff').is_nan()
        ).with_column(
            pl.when(pl.col('p_val') == 0).then(300).otherwise(-pl.col('p_val').log10()).alias('log_p_val')
        ).select([
            'log_p_val',
            pl.col('log_p_val').pow(2).alias('log_p_val2'),
            named_conditions[group_comp[1]][1].alias('subgroup'),
            (pl.col('coeff')*pl.col(f'{ethnicity}_coeff') > 0).cast(int).alias('outcome')
        ]).to_pandas()
        mod = statsmodels.formula.api.logit(
            'outcome ~ log_p_val + log_p_val2 + subgroup', category_df
        ).fit(disp=0)
        p_val = mod.pvalues['subgroup[T.True]']/2
        param = mod.params['subgroup[T.True]']
        if mod.params['subgroup[T.True]'] <= 0:
            p_val += 0.5
        print(f'{ethnicity}: P replication coeff({named_conditions[group_comp[1]][0]}) > 0 in '
              f'{named_conditions[group_comp[0]][0]}: '
              f'{p_val:.3} {"**" if p_val*25 < 0.05 else "*" if p_val < 0.05 else ""} '
              f'param: {param:.3} '
         )

df = df.with_column((-pl.col('p_val').log10()).alias('log_p_val'))

print('plotting supporting figs', flush=True)
figs = []
for ethnicity_num, ethnicity in enumerate(other_ethnicities):
    bin_bounds = [-np.log10(5e-8), 10, 20, 40, 80, np.inf]
    category_names = []
    ys = []
    counts = []
    upper_whiskers = []
    lower_whiskers = []
    for lower_bound, upper_bound in zip(bin_bounds[:-1], bin_bounds[1:]):
        for finemapper, condition in named_conditions:
            name = '['
            if isinstance(lower_bound, float):
                name += f'{lower_bound:.3}'
            else:
                name += str(lower_bound)
            name += ', '
            if isinstance(upper_bound, float):
                name += f'{upper_bound:.3}'
            else:
                name += str(upper_bound)
            name += ')'
            category_names.append((name, finemapper))

            category_df = df.filter(
                (pl.col('log_p_val') >= lower_bound) &
                (pl.col('log_p_val') < upper_bound) &
                condition &
                ~pl.col(f'{ethnicity}_coeff').is_nan()
            )
            bin_size = category_df.shape[0]
            successes = category_df.filter(
                (pl.col('coeff')*pl.col(f'{ethnicity}_coeff') > 0)
            ).shape[0]
            if bin_size == 0:
                ys.append(np.nan)
                counts.append(0)
                upper_whiskers.append(np.nan)
                lower_whiskers.append(np.nan)
                continue
            ys.append(successes/bin_size)
            counts.append(bin_size)
            lower, upper = scipy.stats.binomtest(successes, bin_size).proportion_ci(.5)
            upper_whiskers.append(upper)
            lower_whiskers.append(lower)
    cds = bokeh.models.ColumnDataSource(dict(
        category_names=category_names,
        finemappers=[name_cond[0] for name_cond in named_conditions]*(len(upper_whiskers)//len(named_conditions)),
        ys=ys,
        counts=[str(count) for count in counts],
    ))
    fig = bokeh.plotting.figure(
        title=ethnicity.replace("_", " ").title(),
        y_axis_label = 'Percent loci with shared effect direction',
        x_axis_label = 'Discovery p-value',
        x_range = bokeh.models.FactorRange(*category_names),
        width=1600,
        height=800,
        toolbar_location=None,
        min_border_left=20,
        min_border_right=20,
        output_backend='svg',
        outline_line_color=None,
    )
    fig.y_range.start=0.5
    fig.y_range.end=1.05
    fig.yaxis.bounds=(0.5, 1)
    fig.yaxis.minor_tick_line_color=None
    figs.append(fig)
    fig.xaxis.major_tick_line_color=None
    fig.xaxis.major_label_text_color=None
    fig.xgrid.grid_line_color=None
    fig.xaxis.group_text_font_size = '24px'
    fig.xaxis.subgroup_text_font_size = '24px'
    fig.axis.axis_label_text_font_size = '36px'
    fig.axis.major_label_text_font_size = '30px'
    fig.title.text_font_size = '36px'
    kws = dict(
        x='category_names',
        top='ys',
        width=0.9,
        fill_color = bokeh.transform.factor_cmap('category_names', palette=('#cc6c28', '#e8a723', '#588197', '#a48387'), factors=[name_cond[0] for name_cond in named_conditions], start=1, end=2),
        source=cds,
        alpha=0.7
    )
    if ethnicity_num in (0,3):
        kws['legend_group'] = 'finemappers'

    fig.vbar(
        **kws
    )
    cds = bokeh.models.ColumnDataSource(dict(
        x = category_names,
        x_offset = [-15]*len(category_names),
        y = np.array(ys) + 0.005,
        text = [(f'{count:,.0f}' if count < 1000 else human_format(count)) for count in counts]
    ))
    fig.add_layout(bokeh.models.LabelSet(
        x='x',
        x_offset='x_offset',
        y='y',
        text='text',
        source=cds,
        text_font_size='20px',
        text_color='black'
    ))
    if ethnicity_num in (0,3):
        fig.legend.location = 'top_left'
        fig.legend.label_text_font_size='24px'
        fig.legend[0].items.insert(0, fig.legend[0].items[3])
        del fig.legend[0].items[4]

for fig_n in 1,2,4:
    figs[fig_n].yaxis.axis_label = None
    figs[fig_n].yaxis.major_label_text_font_size = '0px'

bokeh.io.export_png(bokeh.layouts.row(figs[:3]), filename=f'{ukb}/post_finemapping/results/replication_by_pval_non_white.png')
bokeh.io.export_svg(bokeh.layouts.row(figs[:3]), filename=f'{ukb}/post_finemapping/results/replication_by_pval_non_white.svg')
bokeh.io.export_png(bokeh.layouts.row(figs[3:]), filename=f'{ukb}/post_finemapping/results/replication_by_pval_white.png')
bokeh.io.export_svg(bokeh.layouts.row(figs[3:]), filename=f'{ukb}/post_finemapping/results/replication_by_pval_white.svg')
