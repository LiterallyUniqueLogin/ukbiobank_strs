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

print('loading data and running calculations', flush=True)

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

def barplot_fig(cats, data, p_groups):
    fig = bokeh.plotting.figure(
        y_axis_label = 'Percent loci with shared effect direction',
        x_range = bokeh.models.FactorRange(*cats),
        width=len(data)//3*200,
        height=800,
        toolbar_location=None,
        y_range=[.5,1.05],
        output_backend='svg'
    )
    fig.xaxis.major_label_text_color = None
    fig.xaxis.major_tick_line_color = None
    fig.yaxis.bounds = [.5, 1]
    arr_data = np.array(data)
    arr_data[arr_data == 0] = np.nan
    condition_names = np.unique(list(cat[1] for cat in cats))
    cds = bokeh.models.ColumnDataSource(dict(
        x=cats,
        top=arr_data,
        subcats=[cat[1] for cat in cats]
    ))
    fig.vbar(
        x='x',
        top='top',
        source=cds,
        width=0.9,
        line_color='grey',
        fill_color = bokeh.transform.factor_cmap('x', palette=bokeh.palettes.Colorblind[4], factors=condition_names, start=1, end=2),
        legend_group = 'subcats',
    )
    fig.legend.location = 'top_left'
    fig.legend[0].items.insert(0, fig.legend[0].items[3])
    del fig.legend[0].items[4]
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.x_range.range_padding = 0.1
    fig.xaxis.major_label_orientation = 1
    fig.xgrid.grid_line_color = None
    y_range = max(data) - 0.5 # since plots are relative to zero
    n_tests = len(p_groups)*5
    '''
    for group in range(len(p_groups)):
        compare_bars(fig, y_range, cats[group*4], cats[group*4+1], data[group*4], data[group*4+1], data[group*4:group*4+2], p_groups[group][0], 1, -1, -1, n_tests)
        compare_bars(fig, y_range, cats[group*4], cats[group*4+2], data[group*4], data[group*4+2], data[group*4:group*4+3], p_groups[group][1], -1, -1, -2 if p_groups[group][0] <= 0.05 else -1, n_tests)
        compare_bars(fig, y_range, cats[group*4], cats[group*4+3], data[group*4], data[group*4+3], data[group*4:group*4+4], p_groups[group][2], 0, 1, sum(p_groups[group][i] <= 0.05 for i in (2,3,4)), n_tests)
        compare_bars(fig, y_range, cats[group*4+1], cats[group*4+3], data[group*4+1], data[group*4+3], data[group*4+1:group*4+4], p_groups[group][3], 0, 0, sum(p_groups[group][i] <= 0.05 for i in (3,4)), n_tests)
        compare_bars(fig, y_range, cats[group*4+2], cats[group*4+3], data[group*4+2], data[group*4+3], data[group*4+2:group*4+4], p_groups[group][4], 0, -1, 1, n_tests)
    '''
    return fig

def compare_bars(fig, y_range, x1, x2, y1, y2, ys, p_val, x1step, x2step, ystep, n_tests):
    if p_val > 0.05 or np.isnan(p_val):
        return
    x1 = [*x1, .15*x1step]
    x2 = [*x2, .15*x2step]
    y_step_size = 0.05*y_range
    if ystep > 0:
        top = max(ys) + ystep*y_step_size
        text_y = top + 0.01*y_step_size
    else:
        top = min(ys) + ystep*y_step_size
        text_y = top - 0.75*y_step_size
    fig.line(x=[x1, x1], y=[y1, top], color='black')
    fig.line(x=[x2, x2], y=[y2, top], color='black')
    fig.line(x=[x1, x2], y=[top, top], color='black')
    text=f'(p={p_val:.1g})'
    '''
    if p_val < 0.05/n_tests:
        text = '** ' + text
    elif p_val < 0.05:
        text = '* ' + text
    '''
    if p_val < 0.05:
        text = '* ' + text

    cds = bokeh.models.ColumnDataSource(dict(
        x=[x1[:2]],
        y=[text_y],
        #y=[top+0.01*y_step_size],
        text=[text]
    ))
    fig.add_layout(bokeh.models.LabelSet(x='x', y='y', x_offset=.05, text='text', source=cds, text_font_size='12px', text_color='black'))

named_conditions = [
    ('genome-wide significant STRs', pl.col('p_val') <= 5e-8),
    ('FINEMAP STRs with PIP >= 0.8', pl.col('finemapped_finemap')),
    ('SuSiE STRs with max alpha >= 0.8', pl.col('finemapped_susie')),
    ('confidently fine-mapped STRs', pl.col('is_causal_STR_candidate'))
]
'''
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
'''

df = df.with_column((-pl.col('p_val').log10()).alias('log_p_val'))

print('plotting supporting figs', flush=True)
# plot replication across fine-mappers stratified by signal power
figs = []
#finemappers = ['only gwsig', 'finemap', 'susie', 'both', 'confident']
#conditions = [pl.col('p_val') <= 5e-8, pl.col('finemapped_finemap'), pl.col('finemapped_susie'), pl.col('finemapped_finemap') & pl.col('finemapped_susie'), pl.col('is_causal_STR_candidate')]
#finemappers = ['only gwsig', 'finemap', 'susie', 'confident']
#conditions = [pl.col('p_val') <= 5e-8, pl.col('finemapped_finemap'), pl.col('finemapped_susie'), pl.col('is_causal_STR_candidate')]
for ethnicity_num, ethnicity in enumerate(other_ethnicities):
    '''
    quantiles = df.filter(pl.col('log_p_val') >= 10).select([
        pl.col('log_p_val').quantile(.25).alias('q1'),
        pl.col('log_p_val').quantile(.5).alias('q2'),
        pl.col('log_p_val').quantile(.75).alias('q3')
    ])
    print(quantiles)
    bin_bounds = [-np.log10(5e-8), 10, quantiles['q1'][0], quantiles['q2'][0], quantiles['q3'][0], np.inf]
    '''
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
        #title=f'Effect direction sharing in the {ethnicity.replace("_", " ").capitalize()} population by p_val and finemapper',
        title=ethnicity.replace("_", " ").title(),
        y_axis_label = 'Percent loci with shared effect direction',
        x_axis_label = 'Discovery p-value',
        x_range = bokeh.models.FactorRange(*category_names),
        #y_range = [0.5, 1],
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
