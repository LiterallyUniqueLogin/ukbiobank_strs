#!/usr/bin/env python3

import argparse

import bokeh.io
import bokeh.models
import bokeh.plotting
import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outdir')
parser.add_argument('table')
args = parser.parse_args()

finemap_pips = pl.scan_csv(
    args.table, sep='\t'
).filter(
    pl.col('p_val') <= 5e-8
).sort(
    'finemap_pip'
).collect()['finemap_pip'].to_numpy()
cdf = np.cumsum(finemap_pips)/np.sum(finemap_pips)
finemap_frac_15_or_less = np.sum(finemap_pips[finemap_pips <= .15])/np.sum(finemap_pips)

fig = bokeh.plotting.figure(
    title='Total contribution of all variants to FINEMAP PIPs, ordered by those PIPs',
    x_range=(0,1),
    y_range=(0,1),
    x_axis_label='vaiant PIP',
    y_axis_label='total contribution by variants with <= PIP',
    width=1200,
    height=1200
)
fig.step(x=finemap_pips, y=cdf, mode='center')
fig.title.text_font_size = '30px'
fig.axis.axis_label_text_font_size = '26px'
fig.axis.major_label_text_font_size = '20px'
#bokeh.io.export_png(fig, filename=f'{args.outdir}/finemap_pip_cdf.png')
#bokeh.io.export_svg(fig, filename=f'{args.outdir}/finemap_pip_cdf.svg')

susie_alphas = pl.scan_csv(
    args.table, sep='\t'
).filter(
    (pl.col('p_val') <= 5e-8) &
    (pl.col('susie_cs') >= 0)
).sort(
    'susie_alpha'
).collect()['susie_alpha'].to_numpy()
cdf = np.cumsum(susie_alphas)/np.sum(susie_alphas)
susie_frac_15_or_less = np.sum(susie_alphas[susie_alphas <= .15])/np.sum(susie_alphas)

fig = bokeh.plotting.figure(
    title='Total contribution of all variants to FINEMAP PIPs, ordered by those PIPs',
    x_range=(0,1),
    y_range=(0,1),
    x_axis_label='vaiant PIP',
    y_axis_label='total contribution by variants with <= PIP',
    width=1200,
    height=1200
)
fig.step(x=susie_alphas, y=cdf, mode='center')
fig.title.text_font_size = '30px'
fig.axis.axis_label_text_font_size = '26px'
fig.axis.major_label_text_font_size = '20px'
#bokeh.io.export_png(fig, filename=f'{args.outdir}/susie_alpha_cdf.png')
#bokeh.io.export_svg(fig, filename=f'{args.outdir}/susie_alpha_cdf.svg')

with open(f'{args.outdir}/snp_vs_str_stats.txt', 'w') as out:
    out.write(
        f'Among vars with p-val <= 5e-8, vars with PIP <= 0.15 make {susie_frac_15_or_less*100:.2f}% '
        f'of all SuSiE PIP contributions and {finemap_frac_15_or_less*100:.2f}% '
        'of all FINEMAP PIP contributions\n'
    )
    out.flush()

    finemap_gwas_sig_by_var_type = pl.scan_csv(
        args.table, sep='\t'
    ).filter(
        pl.col('p_val') <= 5e-8
    ).groupby(
        'is_STR'
    ).agg([
        pl.col('finemap_pip').sum()
    ]).collect()
    finemap_gwas_sig_str_fraction = (
        finemap_gwas_sig_by_var_type.filter('is_STR')['finemap_pip'].to_numpy()[0] /
        finemap_gwas_sig_by_var_type.select(pl.col('finemap_pip').sum())['finemap_pip'].to_numpy()[0]
    )

    susie_gwas_sig_by_var_type = pl.scan_csv(
        args.table, sep='\t'
    ).filter(
        (pl.col('p_val') <= 5e-8) &
        (pl.col('susie_cs') >= 0)
    ).groupby(
        'is_STR'
    ).agg([
        pl.col('susie_alpha').sum()
    ]).collect()
    susie_gwas_sig_str_fraction = (
        susie_gwas_sig_by_var_type.filter('is_STR')['susie_alpha'].to_numpy()[0] /
        susie_gwas_sig_by_var_type.select(pl.col('susie_alpha').sum())['susie_alpha'].to_numpy()[0]
    )

    out.write(
        f'Among vars with p-val <= 5e-8, STRs make {susie_gwas_sig_str_fraction*100:.2f}% '
        f'of all SuSiE PIP contributions and {finemap_gwas_sig_str_fraction*100:.2f}% '
        'of all FINEMAP PIP contributions\n'
    )

    n_susie_credible_sets = pl.scan_csv(
        args.table, sep='\t'
    ).filter(
        (pl.col('susie_cs') >= 0) &
        (pl.col('p_val') <= 5e-8)
    ).distinct(subset=[
        'phenotype',
        'region',
        'susie_cs'
    ]).collect().shape[0]
    n_regions = pl.scan_csv(
        args.table, sep='\t'
    ).distinct(subset=[
        'phenotype',
        'region',
    ]).collect().shape[0]

    out.write(
        f'SuSiE identifies {n_susie_credible_sets} credible sets across all trait-regions, '
        f'averaging to {n_susie_credible_sets/n_regions:.2f} per region\n'
    )

    finemap_confident_by_var_type = pl.scan_csv(
        args.table, sep='\t'
    ).filter(
        (pl.col('p_val') <= 5e-8) &
        (pl.col('finemap_pip') >= 0.8)
    ).groupby(
        'is_STR'
    ).agg([
        pl.count()
    ]).collect()
    finemap_confident_vars = finemap_confident_by_var_type.select(pl.col('count').sum())['count'].to_numpy()[0]
    finemap_confident_str_fraction = (
        finemap_confident_by_var_type.filter('is_STR')['count'].to_numpy()[0] /
        finemap_confident_vars
    )

    susie_confident_by_var_type = pl.scan_csv(
        args.table, sep='\t'
    ).filter(
        (pl.col('p_val') <= 5e-8) &
        (pl.col('susie_alpha') >= 0.8) &
        (pl.col('susie_cs') >= 0)
    ).groupby(
        'is_STR'
    ).agg([
        pl.count()
    ]).collect()
    susie_confident_vars = susie_confident_by_var_type.select(pl.col('count').sum())['count'].to_numpy()[0]
    susie_confident_str_fraction = (
        susie_confident_by_var_type.filter('is_STR')['count'].to_numpy()[0] /
        susie_confident_vars
    )

    out.write(
        f'Among vars with p-val <= 5e-8, STRs make {susie_confident_str_fraction*100:.2f}% '
        f'of all SuSiE variants with PIP >= 0.8 (out of {susie_confident_vars}) '
        f'and {finemap_confident_str_fraction*100:.2f}% '
        f'of all FINEMAP variants with PIP >= 0.8 (out of {finemap_confident_vars})\n'
    )
    out.flush()

    both_confident_by_var_type = pl.scan_csv(
        args.table, sep='\t'
    ).filter(
        (pl.col('p_val') <= 5e-8) &
        (pl.col('susie_alpha') >= 0.8) &
        (pl.col('susie_cs') >= 0) &
        (pl.col('finemap_pip') >= 0.8)
    ).groupby(
        'is_STR'
    ).agg([
        pl.count()
    ]).collect()
    both_confident_vars = both_confident_by_var_type.select(pl.col('count').sum())['count'].to_numpy()[0]
    both_confident_str_fraction = (
        both_confident_by_var_type.filter('is_STR')['count'].to_numpy()[0] /
        both_confident_vars
    )

    out.write(
        f'Among vars with p-val <= 5e-8, STRs make {both_confident_str_fraction*100:.2f}% '
        f'of all variants with both SuSiE and FINEMAP PIP >= 0.8 (out of {both_confident_vars})\n'
    )

    finemap_by_gwas_sig = pl.scan_csv(
        args.table, sep='\t'
    ).with_column(
        (pl.col('p_val') <= 5e-8).alias('gwas_sig')
    ).groupby(
        'gwas_sig'
    ).agg([
        pl.col('finemap_pip').sum()
    ]).collect()
    finemap_gwas_sig_fraction = (
        finemap_by_gwas_sig.filter('gwas_sig')['finemap_pip'].to_numpy()[0] /
        finemap_by_gwas_sig.select(pl.col('finemap_pip').sum())['finemap_pip'].to_numpy()[0]
    )

    susie_by_gwas_sig = pl.scan_csv(
        args.table, sep='\t'
    ).filter(
        (pl.col('susie_cs') >= 0)
    ).with_column(
        (pl.col('p_val') <= 5e-8).alias('gwas_sig')
    ).groupby(
        'gwas_sig'
    ).agg([
        pl.col('susie_alpha').sum()
    ]).collect()
    susie_gwas_sig_fraction = (
        susie_by_gwas_sig.filter('gwas_sig')['susie_alpha'].to_numpy()[0] /
        susie_by_gwas_sig.select(pl.col('susie_alpha').sum())['susie_alpha'].to_numpy()[0]
    )

    out.write(
        f'Of note, vars with p-val <= 5e-8 make {susie_gwas_sig_fraction*100:.2f}% '
        f'of all SuSiE PIP contributions and {finemap_gwas_sig_fraction*100:.2f}% '
        'of all FINEMAP PIP contributions\n'
    )
