#!/usr/bin/env python3

import argparse
import glob
import os
import os.path
import pathlib
import time

import bokeh.models
import bokeh.models.tickers
import bokeh.io
import bokeh.layouts
import bokeh.plotting
import bokeh.util.hex
import matplotlib.cm
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn
import scipy.stats.kde
import upsetplot

import graphing_utils
import phenotypes

ukb = os.environ['UKB']
other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']

corr_cutoff = .8
p_val_thresh = 5e-8

def all_regions(phenotype):
    return zip(*list(pl.read_csv(
        f'{ukb}/signals/regions/{phenotype}.tab',
        sep='\t'
    ).select([
        (pl.col('chrom').cast(str) + '_' + pl.col('start').cast(str) + '_' + pl.col('end')).alias('region'),
        pl.col('chrom')
    ]).distinct().filter(~(
        (
            (pl.lit(phenotype) == 'total_bilirubin') &
            (pl.col('region') == '12_19976272_22524428')
        ) |
        (
            (pl.lit(phenotype) == 'alkaline_phosphatase') &
            (pl.col('region') == '1_19430673_24309348')
        ) |
        (
            (pl.lit(phenotype) == 'urate') &
            (pl.col('region') == '4_19976272_22524428')
        )
    )).to_dict(False).values()))

def load_susie(results_regions_dir, colnames_regions_dir = None, regions = None, phenotype=None, original=False, return_corrs = False, only_L = None):
    assert (phenotype is not None) + (regions is not None) == 1
    assert (colnames_regions_dir is not None) + original <= 1
    dfs = []
    unconverged_regions = []
    underexplored_regions = []
    unfinished_regions = []
    only_L_skips = []
    min_abs_corrs = []

    if regions is None:
        regions = all_regions(phenotype)
    for (region, chrom) in regions:
        if phenotype:
            print(f'Loading susie {phenotype} region {region}', flush=True)
        if colnames_regions_dir:
            if not os.path.exists(f'{colnames_regions_dir}/{region}'):
                print(f'{colnames_regions_dir}/{region}')
                assert False
        if os.path.exists(f'{results_regions_dir}/{region}/no_strs'):
            continue

        converged_fname = f'{results_regions_dir}/{region}/converged.txt'
        if not os.path.exists(converged_fname):
            unfinished_regions.append(region)
            continue
        with open(converged_fname) as converged_file:
            if not next(converged_file).strip() == 'TRUE':
                unconverged_regions.append(region)
                continue

        if original:
            colnames_fname = f'{results_regions_dir}/{region}/colnames.txt.normal_run'
            if not os.path.exists(colnames_fname):
                colnames_fname = f'{results_regions_dir}/{region}/colnames.txt'
        elif not colnames_regions_dir:
            colnames_fname = f'{results_regions_dir}/{region}/colnames.txt'
        else:
            colnames_fname = f'{colnames_regions_dir}/{region}/colnames.txt'
        if not os.path.exists(colnames_fname) and not original:
            colnames_fname = f'{colnames_fname}.normal_run'
        if not os.path.exists(colnames_fname):
            print(colnames_fname)
            assert False
        with open(colnames_fname) as var_file:
            susie_vars = [line.strip() for line in var_file if line.strip()]

        alphas = np.genfromtxt(
            f'{results_regions_dir}/{region}/alpha.tab',
            delimiter='\t',
        ).T
        if only_L is not None and only_L != alphas.shape[1]:
            only_L_skips.append(region)
            continue
        estimated_signal_vars = np.genfromtxt(
            f'{results_regions_dir}/{region}/V.tab',
            delimiter='\t'
        )

        n_alphas = alphas.shape[1]
        susie_pips=1-np.prod(1-alphas[:, estimated_signal_vars >= 1e-9], axis=1)
        if not susie_pips.shape[0] == len(susie_vars):
            print(results_regions_dir, colnames_fname)
            assert False
        susie_idx = np.arange(len(susie_vars)) + 1
        susie_df = pl.DataFrame({
            'varname': susie_vars,
            'susie_pip': susie_pips,
            'susie_alpha': np.zeros(len(susie_vars)),
            'susie_cs': [-1]*len(susie_vars),
            'susie_idx': susie_idx,
            **{ f'alpha_{i}': alphas[:, i] for i in range(n_alphas) }
        }).lazy()

        df = susie_df.with_columns([
            pl.lit(region).alias('region'),
            pl.lit(chrom).alias('chrom').cast(int),
            pl.col('susie_pip').cast(float).alias('susie_pip')
        ]).sort('susie_idx')

        real_cs_count = 0
        for cs_fname in glob.glob(f'{results_regions_dir}/{region}/cs*.txt'):
            cs_id = int(cs_fname.split('cs')[-1].split('.')[0])
            with open(cs_fname) as cs_file:
                # susie uses 1 based indexing, python uses 0
                # make sure cs idxs are in increasing order
                cs_susie_idx = np.array([int(idx) for idx in next(cs_file).strip().split()])
                assert np.all(cs_susie_idx[1:] - cs_susie_idx[:-1] > 0)
                cs_susie_idx = pl.Series('cs_susie_idx', cs_susie_idx)
                next(cs_file) # skip cs credibility
                min_abs_corr, _, _ = [float(idx) for idx in next(cs_file).strip().split()]
            min_abs_corrs.append(min_abs_corr)
            df = df.with_column(
                pl.when(pl.col('susie_idx').is_in(cs_susie_idx))
                  .then(pl.when(pl.col(f'alpha_{cs_id-1}') > pl.col('susie_alpha'))
                          .then(pl.col(f'alpha_{cs_id-1}'))
                          .otherwise(pl.col('susie_alpha')))
                  .otherwise(pl.col('susie_alpha'))
                  .alias('susie_alpha')
            )
            if min_abs_corr < corr_cutoff:
                continue
            if estimated_signal_vars[cs_id-1] <= 1e-9:
                print(f'CS {cs_id} in region {region} has a pure CS with negligible signal, exiting')
                assert False
            real_cs_count += 1
            if real_cs_count == 50:
                underexplored_regions.append(region)
            # could worry about variants being in multiple CSes
            df = df.with_column(
                pl.when(pl.col('susie_idx').is_in(cs_susie_idx))
                  .then(cs_id)
                  .otherwise(pl.col('susie_cs'))
                  .alias('susie_cs')
            )
        dfs.append(df)

    print('unconverged_regions: ', unconverged_regions)
    print('underexplored_regions: ', underexplored_regions)
    print('unfinished_regions: ', unfinished_regions)
    if only_L is not None:
        print('only_L skips:', only_L_skips)

    dfs = [df.select(
        pl.col('*').exclude('^alpha.*$')
    ) for df in dfs]

    out_df = pl.concat(dfs).drop(['susie_idx'])
    if not return_corrs:
        return out_df
    else:
        return (out_df, min_abs_corrs)

def load_finemap(results_regions_dir, regions = None, phenotype=None):
    assert (phenotype is not None) + (regions is not None) == 1
    unfinished_regions = []
    underexplored_regions = []
    dfs = []

    if regions is None:
        regions = all_regions(phenotype)
    for (region, chrom) in regions:
        if phenotype:
            print(f'Loading finemap {phenotype} region {region}', flush=True)
        if os.path.exists(f'{results_regions_dir}/{region}/no_strs'):
            continue
        if not os.path.exists(f'{results_regions_dir}/{region}/finemap_output.snp') or os.stat(f'{results_regions_dir}/{region}/finemap_output.snp').st_size == 0:
            unfinished_regions.append(region)
            continue
        with open(f'{results_regions_dir}/{region}/finemap_output.log_sss') as log:
            found_n_causal = False
            for line in log:
                if 'n-causal' not in line:
                    continue
                found_n_causal = True
                n_causal = int(line.split()[-1])
                if os.path.exists(f'{results_regions_dir}/{region}/finemap_output.cred{n_causal}'):
                    underexplored_regions.append(region)
                break
            assert found_n_causal

        df = pl.scan_csv(
            f'{results_regions_dir}/{region}/finemap_output.snp',
            sep=' '
        ).select([
            pl.col('rsid').alias('varname'),
            pl.col('prob').alias('finemap_pip'),
            pl.lit(region).alias('region'),
            pl.lit(chrom).alias('chrom').cast(int)
        ])
        dfs.append(df)

    print('unfinished_regions: ', unfinished_regions)
    print('underexplored_regions: ', underexplored_regions)
    return pl.concat(dfs)


def mpv_comparison():
    fname = f'{ukb}/association/results/mean_platelet_volume/my_str/results.tab'
    with open(fname) as tsv:
        header = tsv.readline().strip()
    str_assocs = pl.scan_csv(
        fname,
        sep='\t',
        skip_rows=1,
        has_header=False,
        with_column_names = lambda _: header.replace('0.05_significance_CI', 'foo', 1).replace('5e-8_significance_CI', 'bar', 1).split('\t') # these duplicate column names won't be used anyway
    ).select([
        'chrom',
        'pos',
        ('STR_' + pl.col('pos').cast(str)).alias('varname'),
        pl.lit(True).alias('is_STR'),
        pl.col(f'p_mean_platelet_volume').alias('p_val'),
    ])

    snp_assocs = pl.scan_csv(
        f'{ukb}/association/results/mean_platelet_volume/plink_snp/results.tab',
        sep='\t',
        null_values='NA',
    ).select([
        pl.col('#CHROM').alias('chrom'),
        pl.col('POS').alias('pos'),
        ('SNP_' + pl.col('POS').cast(str) + '_' + pl.col('REF') + '_' +  pl.col('ALT')).alias('varname'),
        pl.lit(False).alias('is_STR'),
        pl.col('P').alias('p_val'),
    ])
    assocs = pl.concat([snp_assocs, str_assocs])

    print('original SuSiE')
    original_susies = load_susie(
        f'{ukb}/finemapping/susie_results/mean_platelet_volume',
        phenotype='mean_platelet_volume',
        original=True,
        only_L=10
    )
    print('tol SuSiE')
    tol_susies = load_susie(
        f'{ukb}/finemapping/susie_results/mean_platelet_volume_tol_0.0001',
        phenotype='mean_platelet_volume',
        colnames_regions_dir=f'{ukb}/finemapping/susie_results/mean_platelet_volume'
    )
    print('res_var 0.8 SuSiE')
    res_var_susies_8 = load_susie(
        f'{ukb}/finemapping/susie_results/mean_platelet_volume_res_var_0.8',
        phenotype='mean_platelet_volume',
        colnames_regions_dir=f'{ukb}/finemapping/susie_results/mean_platelet_volume'
    )
    print('res_var 0.95 SuSiE')
    res_var_susies_95 = load_susie(
        f'{ukb}/finemapping/susie_results/mean_platelet_volume_res_var_0.95',
        phenotype='mean_platelet_volume',
        colnames_regions_dir=f'{ukb}/finemapping/susie_results/mean_platelet_volume'
    )
    print('prior_var 0.2 SuSiE')
    prior_var_susies_02 = load_susie(
        f'{ukb}/finemapping/susie_results/mean_platelet_volume_prior_var_0.2',
        phenotype='mean_platelet_volume',
        colnames_regions_dir=f'{ukb}/finemapping/susie_results/mean_platelet_volume'
    )
    print('prior_var 0.0005 SuSiE')
    prior_var_susies_0005 = load_susie(
        f'{ukb}/finemapping/susie_results/mean_platelet_volume_prior_var_0.0005',
        phenotype='mean_platelet_volume',
        colnames_regions_dir=f'{ukb}/finemapping/susie_results/mean_platelet_volume'
    )
    print('original FINEMAP')
    original_finemaps = load_finemap(f'{ukb}/finemapping/finemap_results/mean_platelet_volume', phenotype='mean_platelet_volume')
    print('total prob FINEMAP')
    total_prob_finemaps        = load_finemap(f'{ukb}/finemapping/finemap_results/mean_platelet_volume.total_prob_4', phenotype='mean_platelet_volume')
    print('prior std FINEMAP derived')
    prior_std_derived_finemaps = load_finemap(f'{ukb}/finemapping/finemap_results/mean_platelet_volume.prior_std_0.0224', phenotype='mean_platelet_volume')
    print('conv tol FINEMAP')
    conv_tol_finemaps          = load_finemap(f'{ukb}/finemapping/finemap_results/mean_platelet_volume.prob_conv_sss_tol_0.0001', phenotype='mean_platelet_volume')
    print('mac FINEMAP')
    mac_finemaps               = load_finemap(f'{ukb}/finemapping/finemap_results_mac_100/mean_platelet_volume', phenotype='mean_platelet_volume')
    print('p threshold FINEMAP')
    p_threshold_finemaps      = load_finemap(f'{ukb}/finemapping/finemap_results_threshold_0.0005/mean_platelet_volume', phenotype='mean_platelet_volume')

    print('Collecting ... ', end='', flush=True)
    start = time.time()
    pheno_df = original_finemaps.join(
        total_prob_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_total_prob'
    ).drop('region_total_prob').join(
        prior_std_derived_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_prior_std_derived'
    ).drop('region_prior_std_dervied').join(
        conv_tol_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_conv_tol'
    ).drop('region_conv_tol').join(
        mac_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_mac'
    ).drop('region_conv_tol').join(
        p_threshold_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_p_thresh'
    ).drop('region_p_thresh').join(
        original_susies,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_extra'
    ).drop('region_extra').join(
        tol_susies,
        how='outer',
        on = ['chrom', 'varname'],
        suffix='_tol'
    ).drop('region_tol').join(
        res_var_susies_8,
        how='outer',
        on = ['chrom', 'varname'],
        suffix='_res_var_8'
    ).drop('region_res_var_8').join(
        res_var_susies_95,
        how='outer',
        on = ['chrom', 'varname'],
        suffix='_res_var_95'
    ).drop('region_res_var_95').join(
        prior_var_susies_02,
        how='outer',
        on = ['chrom', 'varname'],
        suffix='_prior_var_02'
    ).drop('region_prior_var_02').join(
        prior_var_susies_0005,
        how='outer',
        on = ['chrom', 'varname'],
        suffix='_prior_var_0005'
    ).drop('region_prior_var_0005').join(
        assocs,
        how='left',
        on=['chrom', 'varname']
    ).collect()
    print(f'done. Time: {(time.time() - start)/60:.2f}m', flush=True)

    total_df = pheno_df.select([
        'chrom',
        'region',
        'pos',
        'is_STR',
        'varname',
        'p_val',
        'susie_cs',
        'susie_alpha',
        'susie_pip',
        'susie_cs_tol',
        'susie_alpha_tol',
        'susie_pip_tol',
        'susie_cs_res_var_8',
        'susie_alpha_res_var_8',
        'susie_pip_res_var_8',
        'susie_cs_res_var_95',
        'susie_alpha_res_var_95',
        'susie_pip_res_var_95',
        'susie_cs_prior_var_02',
        'susie_alpha_prior_var_02',
        'susie_pip_prior_var_02',
        'susie_cs_prior_var_0005',
        'susie_alpha_prior_var_0005',
        'susie_pip_prior_var_0005',
        'finemap_pip',
        'finemap_pip_total_prob',
        'finemap_pip_prior_std_derived',
        'finemap_pip_conv_tol',
        'finemap_pip_mac',
        'finemap_pip_p_thresh',
    ])

    line_width=18
    min_text_size='30px'
    # concordance graphs
    for suffix, y_label in [
        ('', 'L=10'),
        ('_tol', '10x stricter convergence tolerance'),
        ('_res_var_8', 'residual variance prior = 0.8'),
        ('_res_var_95', 'residual variance prior = 0.95'),
        ('_prior_var_0005', 'prior_var_0005'),
    ]:
        for type_ in ('alpha', 'pip'):
            fig_height = 1200
            fig = bokeh.plotting.figure(
                width=fig_height,
                height=fig_height,
                y_axis_label = ('PIP' if type_ == 'pip' else 'CP') + f' under {y_label}',
                x_axis_label = 'PIP' if type_ == 'pip' else 'CP',
                x_range=[0,1],
                y_range=[0,1],
            )
            fig.title.text_font_size = '36px'
            fig.axis.axis_label_text_font_size = '36px'
            fig.axis.major_label_text_font_size = min_text_size

            graph_df = total_df.filter(
                ~pl.col(f'susie_cs{suffix}').is_null() &
                ~pl.col(f'susie_cs_prior_var_02').is_null()
            )

            if type_ == 'alpha':
                graph_df = graph_df.with_columns([
                    pl.when(
                        pl.col(f'susie_cs{suffix}') >= 0
                    ).then(
                        pl.col(f'susie_alpha{suffix}')
                    ).otherwise(
                        pl.lit(0)
                    ).alias(f'susie_alpha{suffix}'),
                    pl.when(
                        pl.col('susie_cs_prior_var_02') >= 0
                    ).then(
                        pl.col('susie_alpha_prior_var_02')
                    ).otherwise(
                        pl.lit(0)
                    ).alias('susie_alpha_prior_var_02'),
                ])
            else:
                assert type_ == 'pip'
                graph_df = graph_df.with_columns([
                    pl.col(f'susie_pip{suffix}').alias(f'susie_alpha{suffix}'),
                    pl.col('susie_pip_prior_var_02').alias('susie_alpha_prior_var_02'),
                ])
            graph_df = graph_df.filter(
                (pl.col('p_val') <= 5e-8) &
                ((pl.col('susie_alpha_prior_var_02') > 0) | (pl.col(f'susie_alpha{suffix}') > 0))
            )

            n_rects = 100
            graph_df = graph_df.select([
                (pl.col(f'susie_alpha{suffix}')*n_rects).floor()/n_rects,
                (pl.col('susie_alpha_prior_var_02')*n_rects).floor()/n_rects
            ]).groupby(['susie_alpha_prior_var_02', f'susie_alpha{suffix}']).agg(pl.count())

            palette = [linear_int_interpolate((134,204,195), (9,41,46), i/254) for i in range(-1, 255)]
            cmap = bokeh.transform.log_cmap(
                'count',
                palette = palette,
                low=1,
                high=max(graph_df['count'].to_numpy()),
                low_color=(255, 255, 255)
            )
            color_mapper = bokeh.models.LogColorMapper(
                palette = palette,
                low=1,
                high=max(graph_df['count'].to_numpy())
            )

            cds = bokeh.models.ColumnDataSource(dict(
                left=graph_df['susie_alpha_prior_var_02'].to_numpy(), right=graph_df['susie_alpha_prior_var_02'].to_numpy() + 1/n_rects,
                bottom=graph_df[f'susie_alpha{suffix}'].to_numpy(), top=graph_df[f'susie_alpha{suffix}'].to_numpy() + 1/n_rects,
                count=graph_df['count'].to_numpy()
            ))
            fig.quad(left='left', right='right', bottom='bottom', top='top', source=cds, fill_color=cmap, line_width=0)

            thresh = .05
            n_both = graph_df.filter(
                (pl.col(f'susie_alpha{suffix}') >= 1-thresh) &
                (pl.col('susie_alpha_prior_var_02') >= 1-thresh)
            )['count'].sum()

            n_only_finemap = graph_df.filter(
                (pl.col(f'susie_alpha{suffix}') >= 1-thresh) &
                (pl.col('susie_alpha_prior_var_02') <= thresh)
            )['count'].sum()

            n_only_susie = graph_df.filter(
                (pl.col('susie_alpha_prior_var_02') >= 1-thresh) &
                (pl.col(f'susie_alpha{suffix}') <= thresh)
            )['count'].sum()

            xs = np.arange(0, 1, 0.0001)
            fig.line(xs, xs, line_width=18)

            fig.line([1-thresh, 1-thresh], [1-thresh, 1], color='orange', width=line_width)
            fig.line([1-thresh, 1], [1-thresh, 1-thresh], color='orange', width=line_width)
            fig.add_layout(bokeh.models.Title(
                text=f'# variants: {n_both}', align='right',
                text_font_size=min_text_size
            ), 'above')

            fig.line([1-thresh, 1], [thresh, thresh], color='orange', width=line_width)
            fig.line([1-thresh, 1-thresh], [0, thresh], color='orange', width=line_width)
            fig.add_layout(bokeh.models.Title(
                text=f'# variants: {n_only_susie}', align='right',
                text_font_size=min_text_size
            ), 'right')

            fig.line([thresh, thresh], [1-thresh, 1], color='orange', width=line_width)
            fig.line([0, thresh], [1-thresh, 1-thresh], color='orange', width=line_width)
            fig.add_layout(bokeh.models.Title(
                text=f'# variants: {n_only_finemap}', align='left',
                text_font_size=min_text_size
            ), 'above')

            color_bar = bokeh.models.ColorBar(
                title='Number of loci',
                color_mapper = color_mapper,
                width=35,
                height=fig_height//3,
                title_text_font_size='30px',
                ticker = bokeh.models.tickers.FixedTicker(ticks=[0, 10, 100, 1000, 10000, 100000]),
                formatter = bokeh.models.formatters.NumeralTickFormatter(format = '0a'),
                location='bottom_right',
                major_label_text_font_size = min_text_size
            )
            fig.add_layout(color_bar, 'right')
            graphing_utils.resize(fig, 5000/1200, legend=False)
            fig.toolbar_location = None
            fig.background_fill_color = None
            fig.border_fill_color = None
            fig.grid.grid_line_color=None
            bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/consistency/mpv_SuSiE_consistency_{type_}{suffix}.png')
            bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/consistency/mpv_SuSiE_consistency_{type_}{suffix}.svg')

    for suffix, y_label in [
        ('conv_tol', '10x stricter convergence tolerance'),
        ('total_prob', 'assumption of 4 causal variants per region'),
        ('prior_std_derived', 'Avg. 0.05% total var explained per causal variant'),
        ('mac', 'non-major-allele dosage >=100 threshold'),
        ('p_thresh', 'p-value <= 5e-4 threshold')
    ]:
        fig = bokeh.plotting.figure(
            width=1200,
            height=1200,
            y_axis_label = f'CP under {y_label}',
            x_axis_label = 'SuSiE CP',
            x_range=[0,1],
            y_range=[0,1],
            output_backend='svg'
        )
        fig.title.text_font_size = '30px'
        fig.axis.axis_label_text_font_size = '30px'
        fig.axis.major_label_text_font_size = min_text_size

        graph_df = total_df.filter(
            (pl.col('p_val') <= 5e-8) &
            ((pl.col('finemap_pip') > 0) | (pl.col(f'finemap_pip_{suffix}') > 0))
        )

        n_rects = 100
        graph_df = graph_df.select([
            (pl.col(f'finemap_pip_{suffix}')*n_rects).floor()/n_rects,
            (pl.col('finemap_pip')*n_rects).floor()/n_rects
        ]).groupby(['finemap_pip', f'finemap_pip_{suffix}']).agg(pl.count())

        palette = [linear_int_interpolate((134,204,195), (9,41,46), i/254) for i in range(-1, 255)]
        cmap = bokeh.transform.log_cmap(
            'count',
            palette = palette,
            low=1,
            high=max(graph_df['count'].to_numpy()),
            low_color=(255, 255, 255)
        )
        color_mapper = bokeh.models.LogColorMapper(
            palette = palette,
            low=1,
            high=max(graph_df['count'].to_numpy())
        )

        cds = bokeh.models.ColumnDataSource(dict(
            left=graph_df['finemap_pip'].to_numpy(), right=graph_df['finemap_pip'].to_numpy() + 1/n_rects,
            bottom=graph_df[f'finemap_pip_{suffix}'].to_numpy(), top=graph_df[f'finemap_pip_{suffix}'].to_numpy() + 1/n_rects,
            count=graph_df['count'].to_numpy()
        ))
        fig.quad(left='left', right='right', bottom='bottom', top='top', source=cds, fill_color=cmap, line_width=0)

        thresh = .05
        n_both = graph_df.filter(
            (pl.col(f'finemap_pip_{suffix}') >= 1-thresh) &
            (pl.col('finemap_pip') >= 1-thresh)
        )['count'].sum()

        n_only_finemap = graph_df.filter(
            (pl.col(f'finemap_pip_{suffix}') >= 1-thresh) &
            (pl.col('finemap_pip') <= thresh)
        )['count'].sum()

        n_only_susie = graph_df.filter(
            (pl.col('finemap_pip') >= 1-thresh) &
            (pl.col(f'finemap_pip_{suffix}') <= thresh)
        )['count'].sum()

        xs = np.arange(0, 1, 0.0001)
        fig.line(xs, xs, line_width=line_width)

        fig.line([1-thresh, 1-thresh], [1-thresh, 1], color='orange', width=line_width)
        fig.line([1-thresh, 1], [1-thresh, 1-thresh], color='orange', width=line_width)
        fig.add_layout(bokeh.models.Title(
            text=f'# variants: {n_both}', align='right',
            text_font_size=min_text_size
        ), 'above')

        fig.line([1-thresh, 1], [thresh, thresh], color='orange', width=line_width)
        fig.line([1-thresh, 1-thresh], [0, thresh], color='orange', width=line_width)
        fig.add_layout(bokeh.models.Title(
            text=f'# variants: {n_only_susie}', align='right',
            text_font_size=min_text_size
        ), 'right')

        fig.line([thresh, thresh], [1-thresh, 1], color='orange', width=line_width)
        fig.line([0, thresh], [1-thresh, 1-thresh], color='orange', width=line_width)
        fig.add_layout(bokeh.models.Title(
            text=f'# variants: {n_only_finemap}', align='left',
            text_font_size=min_text_size
        ), 'above')

        color_bar = bokeh.models.ColorBar(color_mapper = color_mapper, width=70, major_label_text_font_size = min_text_size)
        fig.add_layout(color_bar, 'right')
        graphing_utils.resize(fig, 5000/1200, legend=False)
        fig.toolbar_location = None
        fig.background_fill_color = None
        fig.border_fill_color = None
        fig.grid.grid_line_color=None
        bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/consistency/mpv_FINEMAP_consistency_{suffix}.png')
        bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/consistency/mpv_FINEMAP_consistency_{suffix}.svg')




        '''
                thresh = .05
        both_df = cs_STRs.filter(
            (pl.col(pip_col) >= 1-thresh) &
            (pl.col(f'{pip_col}_{suffix}') >= 1-thresh)
        )
        n_both = both_df.shape[0]
        avg_both_other_pip = both_df.select(pl.col(other_pip_col).mean())[other_pip_col].to_numpy()[0]

        not_rep_df = cs_STRs.filter(
            (pl.col(pip_col) >= 1-thresh) &
            (pl.col(f'{pip_col}_{suffix}') <= thresh)
        )
        n_not_rep = not_rep_df.shape[0]
        avg_not_rep_other_pip = not_rep_df.select(pl.col(other_pip_col).mean())[other_pip_col].to_numpy()[0]

        new_df = cs_STRs.filter(
            (pl.col(pip_col) <= thresh) &
            (pl.col(f'{pip_col}_{suffix}') >= 1-thresh)
        )
        n_new = new_df.shape[0]
        avg_new_other_pip = new_df.select(pl.col(other_pip_col).mean())[other_pip_col].to_numpy()[0]

        fig = bokeh.plotting.figure(
            width=1200,
            height=1200,
            y_axis_label = f'PIP under {y_label}',
            x_axis_label = 'original PIP',
            x_range=[0,1],
            y_range=[0,1],
            title=f'{mapper} metaparameter comparison'
        )
        fig.title.text_font_size = '30px'
        fig.axis.axis_label_text_font_size = '26px'
        fig.axis.major_label_text_font_size = '20px'

        palette = [
            linear_int_interpolate((111,107,237), (219,46,40), i/254) for i in range(-1, 255)
        ]
        color_mapper = bokeh.models.LinearColorMapper(
            palette = palette,
            low=0,
            high=1
        )
        fig.circle(
            cs_STRs[pip_col].to_numpy(),
            cs_STRs[f'{pip_col}_{suffix}'].to_numpy(),
            size = -np.log10(cs_STRs['p_val'].to_numpy())/7.5,
            alpha = 0.25,
            color=[palette[int(step)] for step in cs_STRs[other_pip_col].to_numpy()*255]
        )

        xs = np.arange(0, 1, 0.0001)
        fig.line(
            xs,
            xs,
            line_dash='dashed'
        )
        fig.quad(
            left=[1-thresh],
            right=[1],
            bottom=[1-thresh],
            top=[1],
            color='orange',
            alpha=0.25
        )
        fig.add_layout(bokeh.models.Title(
            text=f'# STRs: {n_both}, avg {other_label}: {avg_both_other_pip:.2}', align='right',
            text_font_size='18px'
        ), 'above')
        fig.quad(
            left=[1-thresh],
            right=[1],
            bottom=[0],
            top=[thresh],
            color='orange',
            alpha=0.25
        )
        fig.add_layout(bokeh.models.Title(
            text=f'# STRs: {n_not_rep}, avg {other_label}: {avg_not_rep_other_pip:.2}', align='right',
            text_font_size='18px'
        ), 'right')

        fig.quad(
            left=[0],
            right=[thresh],
            bottom=[1-thresh],
            top=[1],
            color='orange',
            alpha=0.25
        )
        fig.add_layout(bokeh.models.Title(
            text=f'# STRs: {n_new}, avg {other_label}: {avg_new_other_pip:.2}', align='left',
            text_font_size='18px'
        ), 'above')

        color_bar = bokeh.models.ColorBar(
            color_mapper = color_mapper,
            width=70,
            title_text_font_size = '26px',
            title=other_label,
            major_label_text_font_size = '20px'
        )
        fig.add_layout(color_bar, 'right')

        fig.toolbar_location = None
        fig.background_fill_color = None
        fig.border_fill_color = None
        fig.grid.grid_line_color=None
        bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/consistency/mpv_{mapper.lower()}_consistency_{suffix}.png')
        bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/consistency/mpv_{mapper.lower()}_consistency_{suffix}.svg')
        '''

def get_putatively_causal_regions(phenotype):
    causal_df = pl.scan_csv(
        f'{ukb}/post_finemapping/results/validated/putatively_causal_STRs.tab',
        sep='\t'
    ).filter(
        pl.col('phenotype') == phenotype
    )
    '''
    pheno_summaries = []
    for phenotype in phenotypes.phenotypes_in_use:
    '''
    pheno_summary = pl.scan_csv(
        f'{ukb}/finemapping/summary/{phenotype}_table.tab',
        sep='\t'
    ).with_column(
        pl.lit(phenotype).alias('phenotype')
    ).select([
        'phenotype', 'chrom', 'start_pos', 'signal_region'
    ])
    '''
        pheno_summaries.append(pheno_summary)
    causal_df = causal_df.join(
        pl.concat(pheno_summaries),
    '''
    causal_df = causal_df.join(
        pheno_summary,
        how='left',
        on=['phenotype', 'chrom', 'start_pos']
    ).select([
        'phenotype', 'signal_region', 'chrom'
    ]).distinct().filter(~(
        (
            (pl.col('phenotype') == 'total_bilirubin') &
            (pl.col('signal_region') == '12_19976272_22524428')
        ) |
        (
            (pl.col('phenotype') == 'mean_platelet_volume') &
            (pl.col('signal_region') == '17_2341352_2710113')
        ) |
        (
            (pl.col('phenotype') == 'alkaline_phosphatase') &
            (pl.col('signal_region') == '1_19430673_24309348')
        )
    )).select([
        'signal_region', 'chrom'
    ]).collect().to_dict(False)

    #phenos, regions, chroms = [
    regions, chroms = [
        causal_df[col] for col in ('signal_region', 'chrom') #('phenotype', 'signal_region', 'chrom')
    ]
    #return list(zip(phenos, regions, chroms))
    return list(zip(regions, chroms))

def putatively_causal_hits_df(phenotype, should_assert):
    regions = get_putatively_causal_regions(phenotype)
    if len(regions) == 0:
        pathlib.Path(f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab').touch()
        pathlib.Path(f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab.empty').touch()
        return

    fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    with open(fname) as tsv:
        header = tsv.readline().strip()
    str_assocs = pl.scan_csv(
        fname,
        sep='\t',
        skip_rows=1,
        has_header=False,
        with_column_names = lambda _: header.replace('0.05_significance_CI', 'foo', 1).replace('5e-8_significance_CI', 'bar', 1).split('\t') # these duplicate column names won't be used anyway
    ).select([
        'chrom',
        'pos',
        ('STR_' + pl.col('pos').cast(str)).alias('varname'),
        pl.lit(True).alias('is_STR'),
        pl.col(f'p_{phenotype}').alias('p_val'),
    ])

    other_ethnicity_assocs = None
    for ethnicity in other_ethnicities:
        fname = f'{ukb}/association/results_finemapped_only/{ethnicity}/{phenotype}/my_str/results.tab'
        with open(fname) as tsv:
            header = tsv.readline().strip()
        one_other_ethnicity = pl.scan_csv(
            fname,
            sep='\t',
            skip_rows=1,
            has_header=False,
            with_column_names = lambda _: header.replace('0.05_significance_CI', 'foo', 1).replace('5e-8_significance_CI', 'bar', 1).split('\t') # these duplicate column names won't be used anyway
        ).select([
            'chrom',
            'pos',
            pl.lit(True).alias('is_STR'),
            pl.col(f'p_{phenotype}').alias(f'{ethnicity}_p_val'),
        ])
        if other_ethnicity_assocs is None:
            other_ethnicity_assocs = one_other_ethnicity
        else:
            other_ethnicity_assocs = other_ethnicity_assocs.join(
                one_other_ethnicity,
                how='inner',
                on=['chrom', 'pos', 'is_STR']
            )

    snp_assocs = pl.scan_csv(
        f'{ukb}/association/results/{phenotype}/plink_snp/results.tab',
        sep='\t',
        null_values='NA',
    ).select([
        pl.col('#CHROM').alias('chrom'),
        pl.col('POS').alias('pos'),
        ('SNP_' + pl.col('POS').cast(str) + '_' + pl.col('REF') + '_' +  pl.col('ALT')).alias('varname'),
        pl.lit(False).alias('is_STR'),
        pl.col('P').alias('p_val'),
    ])
    assocs = pl.concat([snp_assocs, str_assocs])

    print('original SuSiE')
    original_susies = load_susie(f'{ukb}/finemapping/susie_results/{phenotype}', regions=regions)
    print('best guess SuSiE')
    best_guess_susies = load_susie(f'{ukb}/finemapping/susie_hardcall_results/{phenotype}', regions=regions)
    print('ratio SuSiE')
    ratio_susies = load_susie(
        f'{ukb}/finemapping/susie_results/{phenotype}_snp_str_ratio_4',
        colnames_regions_dir=f'{ukb}/finemapping/susie_results/{phenotype}',
        regions=regions
    )

    print('original FINEMAP')
    original_finemaps          = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}', regions=regions)
    print('ratio FINEMAP')
    ratio_finemaps             = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}.snp_str_ratio_4', regions=regions)
    print('total prob FINEMAP')
    total_prob_finemaps        = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}.total_prob_4', regions=regions)
    print('prior std FINEMAP derived')
    prior_std_derived_finemaps = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}.prior_std_0.0224', regions=regions)
    print('prior std FINEMAP low')
    prior_std_low_finemaps     = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}.prior_std_0.005', regions=regions)
    print('conv tol FINEMAP')
    conv_tol_finemaps          = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}.prob_conv_sss_tol_0.0001', regions=regions)
    print('mac FINEMAP')
    mac_finemaps               = load_finemap(f'{ukb}/finemapping/finemap_results_mac_100/{phenotype}', regions=regions)
    print('p threshold FINEMAP')
    p_threshold_finemaps      = load_finemap(f'{ukb}/finemapping/finemap_results_threshold_0.0005/{phenotype}', regions=regions)

    print('Collecting ... ', end='', flush=True)
    start = time.time()
    pheno_df = original_finemaps.join(
        ratio_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_ratio'
    ).drop('region_ratio').join(
        total_prob_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_total_prob'
    ).drop('region_total_prob').join(
        prior_std_low_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_prior_std_low'
    ).drop('region_prior_std_low').join(
        prior_std_derived_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_prior_std_derived'
    ).drop('region_prior_std_derived').join(
        conv_tol_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_conv_tol'
    ).drop('region_conv_tol').join(
        mac_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_mac'
    ).drop('region_mac').join(
        p_threshold_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_p_thresh'
    ).drop('region_p_thresh').join(
        original_susies,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_extra'
    ).drop('region_extra').join(
        best_guess_susies,
        how='outer',
        on = ['chrom', 'varname'],
        suffix='_best_guess'
    ).drop('region_best_guess').join(
        ratio_susies,
        how='outer',
        on = ['chrom', 'varname'],
        suffix='_ratio'
    ).drop('region_ratio').join(
        assocs,
        how='left',
        on=['chrom', 'varname']
    ).join(
        other_ethnicity_assocs,
        how='left',
        on=['chrom', 'pos', 'is_STR']
    ).with_column(
        pl.lit(phenotype).alias('phenotype')
    ).collect()
    print(f'done. Time: {(time.time() - start)/60:.2f}m', flush=True)

    pheno_df = pheno_df.select([
        'phenotype',
        'chrom',
        'region',
        'pos',
        'is_STR',
        'varname',
        'p_val',
        'susie_cs',
        'susie_alpha',
        'susie_pip',
        'susie_cs_best_guess',
        'susie_alpha_best_guess',
        'susie_pip_best_guess',
        'susie_cs_ratio',
        'susie_alpha_ratio',
        'susie_pip_ratio',
        'finemap_pip',
        'finemap_pip_ratio',
        'finemap_pip_total_prob',
        'finemap_pip_prior_std_derived',
        'finemap_pip_prior_std_low',
        'finemap_pip_conv_tol',
        'finemap_pip_mac',
        'finemap_pip_p_thresh',
        *[f'{ethnicity}_p_val' for ethnicity in other_ethnicities],
    ])
    pheno_df.write_csv(f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab', sep='\t')
    if should_assert:
        assert pheno_df.select((pl.col('region') == '').any().alias('region'))['region'].to_numpy()[0] == False
        assert np.all(~np.isnan(pheno_df['p_val'].to_numpy()))
        assert np.all(~np.isnan(pheno_df['finemap_pip'].to_numpy()))
        assert np.all(np.isnan(pheno_df['susie_alpha'].to_numpy()) == np.isnan(pheno_df['susie_alpha_best_guess'].to_numpy()))
        assert np.all(np.isnan(pheno_df['susie_alpha'].to_numpy()) == np.isnan(pheno_df['susie_alpha_ratio'].to_numpy()))
        assert np.all(1 == pheno_df.groupby(['phenotype', 'chrom', 'varname']).agg([pl.count()]).sort('count')['count'].to_numpy())
        for ethnicity in other_ethnicities:
            assert np.all(~np.isnan(pheno_df.filter('is_STR')[f'{ethnicity}_p_val'].to_numpy()))

def linear_int_interpolate(c1, c2, dist):
    c_new = []
    for coord1, coord2 in zip(c1, c2):
        c_new.append(coord1 + round((coord2 - coord1)*dist))
    return c_new

def plot_upset_plots(total_df):
    for upset_thresh in .8, .9:
        susie_cols = total_df.select([
            pl.col('^susie_alpha.*$'),
        ]).columns
        susie_cols.remove('susie_alpha')
        susie_cols.insert(len(susie_cols), 'susie_alpha')

        finemap_cols = total_df.select([
            pl.col('^finemap_pip.*$')
        ]).columns
        finemap_cols.remove('finemap_pip')
        #finemap_cols.remove('finemap_pip_ratio')
        finemap_cols.insert(len(finemap_cols), 'finemap_pip')

        print('Converting to pandas ... ', flush=True)
        for found_in_default in True, False:
            intermediate_df = total_df.filter(
                pl.col('is_STR') &
                (pl.col('p_val') <= 1e-10)
            ).with_columns([
                pl.col('^susie_alpha.*$') >= upset_thresh,
                pl.col('^finemap_pip.*$') >= upset_thresh
            ]).filter(
                # passes thresh in at least one fine-mapping run
                pl.sum([pl.col(col).cast(int) for col in susie_cols + finemap_cols]) > 0
            )
            if found_in_default:
                intermediate_df = intermediate_df.filter(
                    pl.col('susie_alpha') | pl.col('finemap_pip')
                )
            intermediate_df = intermediate_df.to_pandas()

            #all plot
            upset_df = upsetplot.from_indicators([col for col in susie_cols + finemap_cols], data=intermediate_df)
            upsetplot.UpSet(
                upset_df,
                sort_categories_by=None,
                show_counts=True,
            ).plot()
            plt.suptitle('Fine-mapping conditions breakdown')
            plt.savefig(f'{ukb}/post_finemapping/results/upsets/all_{upset_thresh}_found_in_default_{found_in_default}.png')
            
            for name, upset_cols, other_col, other_name in ('susie', susie_cols, 'finemap_pip', 'FINEMAP'), ('finemap', finemap_cols, 'susie_alpha', 'SuSiE'):
                print('Plotting upset ... ', flush=True)
                upset_df = upsetplot.from_indicators(upset_cols, data=intermediate_df)
                print('done.', flush=True)
                upset = upsetplot.UpSet(
                    upset_df,
                    sort_categories_by=None,
                    show_counts=True,
                    intersection_plot_elements=0
                )
                upset.add_stacked_bars(
                    by=other_col,
                    title=f'Count, selected by default {other_name} run or not',
                    elements=10,
                    colors=matplotlib.cm.Pastel1
                )
                upset.plot()
                plt.suptitle(f'{name} fine-mapping conditions breakdown')
                plt.savefig(f'{ukb}/post_finemapping/results/upsets/{name}_{upset_thresh}_found_in_default_{found_in_default}.png')


def putatively_causal_hits_comparison():
    print('loading pre-generated df ...', end='', flush=True)
    cols = pl.read_csv(
        f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_white_blood_cell_count.tab',
        sep='\t',
        n_rows=1
    ).columns
    total_df = pl.concat([
        pl.read_csv(
            f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab',
            sep='\t',
            dtypes={col: (float if 'cs' not in col else int) for col in cols if 'finemap' in col or 'susie' in col or 'p_val' in col}
        ) for phenotype in phenotypes.phenotypes_in_use
        if not os.path.exists(f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab.empty')
    ]).rename({'susie_cs_hardcall' : 'susie_cs_best_guess', 'susie_alpha_hardcall': 'susie_alpha_best_guess', 'susie_pip_hardcall': 'susie_pip_best_guess'})
    print(' done', flush=True)

    total_df = total_df.filter(~pl.col('finemap_pip').is_null() & ~pl.col('susie_alpha').is_null())

    total_df = total_df.with_columns([
        pl.when(pl.col('susie_cs') > 0).then(pl.col('susie_alpha')).otherwise(0).alias('susie_alpha'),
        pl.when(pl.col('susie_cs_ratio') > 0).then(pl.col('susie_alpha_ratio')).otherwise(0).alias('susie_alpha_ratio'),
        pl.when(pl.col('susie_cs_best_guess') > 0).then(pl.col('susie_alpha_best_guess')).otherwise(0).alias('susie_alpha_best_guess'),
    ])
    total_df.filter(
        pl.col('is_STR') &
        (pl.col('p_val') <= 1e-10) &
        (pl.col('susie_alpha') >= .8) &
        (pl.col('finemap_pip') >= .8)
    ).drop('is_STR').write_csv(f'{ukb}/post_finemapping/intermediate_results/original_causal_STR_candidates.tab', sep='\t')

    susie_cols = total_df.select([
        pl.col('^susie_alpha.*$'),
    ]).columns
    finemap_cols = total_df.select([
        pl.col('^finemap_pip.*$')
    ]).columns
    pass_all_threshes = total_df.filter(
        pl.col('is_STR') &
        (pl.col('p_val') <= 1e-10) &
        (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols if 'ratio' not in col and 'prior_std_low' not in col]) == 8)
    ).drop(['susie_cs', 'susie_cs_best_guess', 'susie_cs_ratio', 'is_STR', 'varname'])
    pass_all_threshes.write_csv(f'{ukb}/post_finemapping/intermediate_results/concordant_causal_STR_candidates.tab', sep='\t')
    total_df.filter(
        pl.col('is_STR') &
        (pl.col('p_val') <= 1e-10) &
        (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols]) == 11)
    ).select(['phenotype', 'region', 'chrom', 'pos', 'p_val']).write_csv(f'{ukb}/post_finemapping/intermediate_results/strictly_concordant_causal_STR_candidates.tab', sep='\t')

    # concordance graphs
    for mapper, meas, suffix, y_label in [
        ('SuSiE', 'alpha', 'ratio', '4x prior on SNPs/Indels'),
        ('SuSiE', 'alpha', 'best_guess', 'best guess genotyping'),
        ('SuSiE', 'pip', 'ratio', '4x prior on SNPs/Indels'),
        ('SuSiE', 'pip', 'best_guess', 'best guess genotyping'),
        ('FINEMAP', 'pip', 'ratio', '4x prior on SNPs/Indels'),
        ('FINEMAP', 'pip', 'conv_tol', '10x stricter convergence tolerance'),
        ('FINEMAP', 'pip', 'total_prob', 'assumption of 4 causal variants per region'),
        ('FINEMAP', 'pip', 'prior_std_low', 'Avg. 0.0025% total var explained per causal variant'),
        ('FINEMAP', 'pip', 'prior_std_derived', 'Avg. 0.05% total var explained per causal variant'),
        ('FINEMAP', 'pip', 'mac', 'non-major-allele dosage >=100 threshold'),
        ('FINEMAP', 'pip', 'p_thresh', 'p-value <= 5e-4 threshold')
    ]:
        figs = []
        for var_type, var_cond in (('STR', pl.col('is_STR')), ('SNP/indel', ~pl.col('is_STR'))):
            if var_type != 'STR' and suffix == 'ratio':
                continue
            if mapper == 'SuSiE':
                pip_col = f'susie_{meas}'
                other_label = 'FINEMAP CP'
                other_pip_col = 'finemap_pip'
            else:
                assert mapper == 'FINEMAP'
                pip_col = 'finemap_pip'
                other_label = 'SuSiE CP'
                other_pip_col = 'susie_alpha'

            graph_df = total_df.filter(
                var_cond &
                (pl.col('p_val') <= 1e-10) &
                ((pl.col(pip_col) > 0) | (pl.col(f'{pip_col}_{suffix}') > 0))
            ).sort(pip_col).with_column(
                pl.max([pl.col('p_val'), 1e-300]).alias('p_val')
            )
            thresh = .025
            both_df = graph_df.filter(
                (pl.col(pip_col) >= 1-thresh) &
                (pl.col(f'{pip_col}_{suffix}') >= 1-thresh)
            )
            n_both = both_df.shape[0]
            avg_both_other_pip = both_df.select(pl.col(other_pip_col).mean())[other_pip_col].to_numpy()[0]

            not_rep_df = graph_df.filter(
                (pl.col(pip_col) >= 1-thresh) &
                (pl.col(f'{pip_col}_{suffix}') <= thresh)
            )
            n_not_rep = not_rep_df.shape[0]
            avg_not_rep_other_pip = not_rep_df.select(pl.col(other_pip_col).mean())[other_pip_col].to_numpy()[0]

            new_df = graph_df.filter(
                (pl.col(pip_col) <= thresh) &
                (pl.col(f'{pip_col}_{suffix}') >= 1-thresh)
            )
            n_new = new_df.shape[0]
            avg_new_other_pip = new_df.select(pl.col(other_pip_col).mean())[other_pip_col].to_numpy()[0]

            fig_height = 1200
            fig_kws = dict(
                width=fig_height,
                height=fig_height,
                x_axis_label = 'original CP',
                x_range=[0,1],
                y_range=[0,1],
                match_aspect=True,
                #title=f'{mapper} metaparameter comparison for {var_type}s'
                output_backend='svg'
            )
            if var_type == 'STR':
                fig_kws['y_axis_label'] = f'CP under {y_label}'
            fig = bokeh.plotting.figure(
                **fig_kws
            )
            fig.title.text_font_size = '36px'
            fig.axis.axis_label_text_font_size = '36px'
            min_text_size = '30px'
            fig.axis.major_label_text_font_size = min_text_size

            palette = [
                linear_int_interpolate((111,107,237), (219,46,40), i/254) for i in range(-1, 255)
            ]
            color_mapper = bokeh.models.LinearColorMapper(
                palette = palette,
                low=0,
                high=1
            )
            fig.circle(
                graph_df[pip_col].to_numpy(),
                graph_df[f'{pip_col}_{suffix}'].to_numpy(),
                size = -np.log10(graph_df['p_val'].to_numpy())/7.5,
                alpha = 0.4,
                color=[palette[int(step)] for step in graph_df[other_pip_col].to_numpy()*255]
            )

            xs = np.arange(0, 1, 0.0001)
            fig.line(
                xs,
                xs,
                line_dash='dashed',
                line_width=6,
                color='black',
                alpha=0.6
            )
            box_alpha = 0.4
            fig.quad(
                left=[1-thresh],
                right=[1],
                bottom=[1-thresh],
                top=[1],
                color='orange',
                alpha=box_alpha
            )
            fig.add_layout(bokeh.models.Title(
                text=f'# {var_type}s: {n_both}, avg {other_label}: {avg_both_other_pip:.2}', align='right',
                text_font_size=min_text_size
            ), 'above')
            fig.quad(
                left=[1-thresh],
                right=[1],
                bottom=[0],
                top=[thresh],
                color='orange',
                alpha=box_alpha
            )
            fig.add_layout(bokeh.models.Title(
                text=f'# {var_type}s: {n_not_rep}, avg {other_label}: {avg_not_rep_other_pip:.2}', align='right',
                text_font_size=min_text_size
            ), 'right')

            if n_new > 5:
                fig.quad(
                    left=[0],
                    right=[thresh],
                    bottom=[1-thresh],
                    top=[1],
                    color='orange',
                    alpha=box_alpha
                )
                fig.add_layout(bokeh.models.Title(
                    text=f'# {var_type}s: {n_new}, avg {other_label}: {avg_new_other_pip:.2}', align='left',
                    text_font_size=min_text_size
                ), 'above')
            fig.toolbar_location = None
            fig.background_fill_color = None
            fig.border_fill_color = None
            fig.grid.grid_line_color=None
            figs.append(fig)

        size_scale_fig = bokeh.plotting.figure(
            width=175,
            height=600,
            y_axis_label = '-log10 p-value',
            y_axis_type='log',
            x_axis_label='Size scale',
            output_backend='svg'
        )
        scales=np.array([10, 20, 40, 80, 160, 300])
        size_scale_fig.yaxis.ticker = bokeh.models.tickers.FixedTicker(ticks=scales)
        size_scale_fig.axis.axis_label_text_font_size = '30px'
        size_scale_fig.axis.major_label_text_font_size = '30px'
        size_scale_fig.toolbar_location = None
        size_scale_fig.background_fill_color = None
        size_scale_fig.border_fill_color = None
        size_scale_fig.grid.grid_line_color = None
        size_scale_fig.xaxis.visible = False

        size_scale_fig.circle(
            np.zeros(scales.shape),
            scales,
            size = scales/7.5,
            color='blue',
            alpha=0.4
        )

        color_bar = bokeh.models.ColorBar(
            title=other_label,
            color_mapper = color_mapper,
            width=35,
            height=fig_height//3,
            major_label_text_font_size = '30px',
            title_text_font_size='30px',
            location='bottom_right'
        )
        figs[-1].add_layout(color_bar, 'right')

        total_fig = bokeh.layouts.row([*figs, size_scale_fig])
        bokeh.io.export_png(total_fig, filename=f'{ukb}/post_finemapping/results/consistency/{mapper.lower()}_consistency_{meas}_{suffix}.png')
        bokeh.io.export_svg(total_fig, filename=f'{ukb}/post_finemapping/results/consistency/{mapper.lower()}_consistency_{meas}_{suffix}.svg')

def first_pass_df(phenotype, should_assert):
    #pheno_dfs = []
    min_abs_corrs = []
    #for count, phenotype in [(0, 'aspartate_aminotransferase'), (1, 'total_bilirubin')]:
    #for count, phenotype in [(0, 'mean_platelet_volume')]:
    #for count, phenotype in enumerate(phenotypes.phenotypes_in_use):
    #    print(f"Loading phenotype #{count+1} ({phenotype})", flush=True)

    fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    with open(fname) as tsv:
        header = tsv.readline().strip()
    str_assocs = pl.scan_csv(
        fname,
        sep='\t',
        skip_rows=1,
        has_header=False,
        with_column_names = lambda _: header.replace('0.05_significance_CI', 'foo', 1).replace('5e-8_significance_CI', 'bar', 1).split('\t') # these duplicate column names won't be used anyway
    ).select([
        'chrom',
        'pos',
        ('STR_' + pl.col('pos').cast(str)).alias('varname'),
        pl.lit(True).alias('is_STR'),
        pl.col(f'p_{phenotype}').alias('p_val'),
        pl.col(f'coeff_{phenotype}').alias('coeff'),
        pl.col(f'se_{phenotype}').alias('se')
    ])

    snp_assocs = pl.scan_csv(
        f'{ukb}/association/results/{phenotype}/plink_snp/results.tab',
        sep='\t',
        null_values='NA',
    ).select([
        pl.col('#CHROM').alias('chrom'),
        pl.col('POS').alias('pos'),
        ('SNP_' + pl.col('POS').cast(str) + '_' + pl.col('REF') + '_' +  pl.col('ALT')).alias('varname'),
        pl.lit(False).alias('is_STR'),
        pl.col('P').alias('p_val'),
        pl.col('BETA').alias('coeff'),
        pl.col('SE').alias('se'),
    ])
    assocs = pl.concat([snp_assocs, str_assocs])

    other_ethnicity_assocs = None
    other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']
    for ethnicity in other_ethnicities:
        fname = f'{ukb}/association/results_finemapped_only/{ethnicity}/{phenotype}/my_str/results.tab'
        with open(fname) as tsv:
            header = tsv.readline().strip()
        one_other_ethnicity = pl.scan_csv(
            fname,
            sep='\t',
            skip_rows=1,
            has_header=False,
            with_column_names = lambda _: header.replace('0.05_significance_CI', 'foo', 1).replace('5e-8_significance_CI', 'bar', 1).split('\t') # these duplicate column names won't be used anyway
        ).select([
            'chrom',
            'pos',
            pl.lit(True).alias('is_STR'),
            pl.col(f'p_{phenotype}').alias(f'{ethnicity}_p_val'),
            pl.col(f'coeff_{phenotype}').alias(f'{ethnicity}_coeff'),
            pl.col(f'se_{phenotype}').alias(f'{ethnicity}_se'),
        ])
        if other_ethnicity_assocs is None:
            other_ethnicity_assocs = one_other_ethnicity
        else:
            other_ethnicity_assocs = other_ethnicity_assocs.join(
                one_other_ethnicity,
                how='inner',
                on=['chrom', 'pos', 'is_STR']
            )

    print('original SuSiE')
    original_susies, pheno_min_abs_corrs = load_susie(f'{ukb}/finemapping/susie_results/{phenotype}', phenotype=phenotype, original=True, return_corrs=True)
    min_abs_corrs.extend(pheno_min_abs_corrs)
    print('original FINEMAP')
    original_finemaps = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}', phenotype=phenotype)

    print('Collecting ... ', end='', flush=True)
    start = time.time()
    pheno_df = original_finemaps.join(
        original_susies,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_extra'
    ).drop('region_extra').join(
        assocs,
        how='left',
        on=['chrom', 'varname']
    ).join(
        other_ethnicity_assocs,
        how='left',
        on=['chrom', 'pos', 'is_STR']
    ).with_column(
        pl.lit(phenotype).alias('phenotype')
    ).collect()
    print(f'done. Time: {(time.time() - start)/60:.2f}m', flush=True)

    # choose col order
    total_df = pheno_df.select([
        'phenotype',
        'chrom',
        'region',
        'pos',
        'is_STR',
        'varname',
        'p_val',
        'coeff',
        'se',
        'susie_pip',
        'susie_cs',
        'susie_alpha',
        'finemap_pip',
        *[f'{ethnicity}_p_val' for ethnicity in other_ethnicities],
        *[f'{ethnicity}_coeff' for ethnicity in other_ethnicities],
        *[f'{ethnicity}_se' for ethnicity in other_ethnicities],
    ])
    total_df.write_csv(f'{ukb}/post_finemapping/intermediate_results/finemapping_all_concordance_{phenotype}.tab', sep='\t')
    np.save(f'{ukb}/post_finemapping/intermediate_results/susie_all_min_abs_corrs_{phenotype}.npy', np.array(min_abs_corrs))
    if should_assert:
        assert total_df.select((pl.col('region') == '').any().alias('region'))['region'].to_numpy()[0] == False
        assert np.all(~np.isnan(total_df['p_val'].to_numpy()))
        assert np.all(~np.isnan(total_df['finemap_pip'].to_numpy()))
        assert np.all(1 == total_df.groupby(['phenotype', 'chrom', 'varname']).agg([pl.count()]).sort('count')['count'].to_numpy())

def first_pass_comparison():
    total_df = pl.concat([
        pl.read_csv(
            f'{ukb}/post_finemapping/intermediate_results/finemapping_all_concordance_{phenotype}.tab',
            sep='\t',
            dtypes={
                **{f'{ethnicity}_p_val': float for ethnicity in other_ethnicities},
                **{f'{ethnicity}_coeff': float for ethnicity in other_ethnicities},
                **{f'{ethnicity}_se': float for ethnicity in other_ethnicities}
            }
        ) for phenotype in phenotypes.phenotypes_in_use
    ])
    min_abs_corrs = np.concatenate([
        np.load(f'{ukb}/post_finemapping/intermediate_results/susie_all_min_abs_corrs_{phenotype}.npy')
        for phenotype in phenotypes.phenotypes_in_use
    ])

    # ----- Stats ----

    print(
        'Total SuSiE CSes',
        # should be the same
        total_df.filter(pl.col('susie_cs') >= 0).unique(subset = ['phenotype', 'region', 'susie_cs']).shape[0],
        np.sum(min_abs_corrs >= corr_cutoff)
    )
    strong_susie_vars = (pl.col('p_val') <= 5e-8) & (pl.col('susie_cs') >= 0) & (pl.col('susie_alpha') >= 0.8)
    print(
        'Total GWsig SuSiE vars PIP >= 0.8, fraction of those that are STRs, (min, max fraction over phenotypes)',
        total_df.filter(strong_susie_vars).shape[0],
        total_df.filter(strong_susie_vars & pl.col('is_STR')).shape[0]/total_df.filter(strong_susie_vars).shape[0],
        total_df.filter(
            strong_susie_vars
        ).groupby([
            'phenotype', 'is_STR'
        ]).agg(pl.count()).with_column(
            pl.col('count').sum().over('phenotype').alias('total')
        ).filter('is_STR').select([
            (pl.col('count')/pl.col('total')).min().alias('min'),
            (pl.col('count')/pl.col('total')).max().alias('max')
        ])
    )
    strong_finemap_vars = (pl.col('p_val') <= 5e-8) & (pl.col('finemap_pip') >= 0.8)
    print(
        'Total GWsig FINEMAP vars PIP >= 0.8, fraction of those that are STRs, (min, max fraction over phenotypes)',
        total_df.filter(strong_finemap_vars).shape[0],
        total_df.filter(strong_finemap_vars & pl.col('is_STR')).shape[0]/total_df.filter(strong_finemap_vars).shape[0],
        total_df.filter(
            strong_finemap_vars
        ).groupby([
            'phenotype', 'is_STR'
        ]).agg(pl.count()).with_column(
            pl.col('count').sum().over('phenotype').alias('total')
        ).filter('is_STR').select([
            (pl.col('count')/pl.col('total')).min().alias('min'),
            (pl.col('count')/pl.col('total')).max().alias('max')
        ])
    )
    strong_both_vars = strong_finemap_vars & strong_susie_vars
    print(
        'Total GWsig FINEMAP and SuSiE vars PIP >= 0.8, fraction of those that are STRs, (min, max fraction over phenotypes)',
        total_df.filter(strong_both_vars).shape[0],
        total_df.filter(strong_both_vars & pl.col('is_STR')).shape[0]/total_df.filter(strong_both_vars).shape[0],
        total_df.filter(
            strong_both_vars
        ).groupby([
            'phenotype', 'is_STR'
        ]).agg(pl.count()).with_column(
            pl.col('count').sum().over('phenotype').alias('total')
        ).filter('is_STR').select([
            (pl.col('count')/pl.col('total')).min().alias('min'),
            (pl.col('count')/pl.col('total')).max().alias('max')
        ])
    )
    
    print(
        'Total contribution of STRs in SuSiE (min, max fraction over phenotypes)',
        total_df.filter(
            (pl.col('p_val') <= 5e-8)
            & (pl.col('susie_cs') >= 0)
        ).groupby('is_STR').agg(
            pl.col('susie_alpha').sum()
        ).with_column(
            pl.col('susie_alpha').sum().alias('total')
        ).filter('is_STR').select(
            (pl.col('susie_alpha')/pl.col('total'))
        ),
        total_df.filter(
            (pl.col('p_val') <= 5e-8)
            & (pl.col('susie_cs') >= 0)
        ).groupby(['phenotype', 'is_STR']).agg(
            pl.col('susie_alpha').sum()
        ).with_column(
            pl.col('susie_alpha').sum().over('phenotype').alias('total')
        ).filter('is_STR').select([
            (pl.col('susie_alpha')/pl.col('total')).min().alias('min'),
            (pl.col('susie_alpha')/pl.col('total')).max().alias('max')
        ])
    )
    print(
        'Total contribution of STRs in FINEMAP (min, max fraction over phenotypes)',
        total_df.filter(
            (pl.col('p_val') <= 5e-8)
            & ~pl.col('finemap_pip').is_null()
        ).groupby('is_STR').agg(
            pl.col('finemap_pip').sum()
        ).with_column(
            pl.col('finemap_pip').sum().alias('total')
        ).filter('is_STR').select(
            (pl.col('finemap_pip')/pl.col('total'))
        ),
        total_df.filter(
            (pl.col('p_val') <= 5e-8)
            & ~pl.col('finemap_pip').is_null()
        ).groupby(['phenotype', 'is_STR']).agg(
            pl.col('finemap_pip').sum()
        ).with_column(
            pl.col('finemap_pip').sum().over('phenotype').alias('total')
        ).filter('is_STR').select([
            (pl.col('finemap_pip')/pl.col('total')).min().alias('min'),
            (pl.col('finemap_pip')/pl.col('total')).max().alias('max')
        ])
    )

    print(
        'Total SuSiE contribution of variants with PIP <= 0.1',
        total_df.filter(
            (pl.col('p_val') <= 5e-8)
            & (pl.col('susie_cs') >= 0)
        ).groupby(
            (pl.col('susie_alpha') < 0.1).alias('low_alpha')
        ).agg(pl.col('susie_alpha').sum()).with_column(
            pl.col('susie_alpha').sum().alias('total')
        ).filter('low_alpha').select(
            pl.col('susie_alpha')/pl.col('total')
        )
    )

    print(
        'Total SuSiE contribution of variants with PIP <= 0.1',
        total_df.filter(
            (pl.col('p_val') <= 5e-8)
            & ~pl.col('finemap_pip').is_null()
        ).groupby(
            (pl.col('finemap_pip') < 0.1).alias('low_alpha')
        ).agg(pl.col('finemap_pip').sum()).with_column(
            pl.col('finemap_pip').sum().alias('total')
        ).filter('low_alpha').select(
            pl.col('finemap_pip')/pl.col('total')
        )
    )

    print(
        'Percentage of SuSiE CSes with FINEMAP PIP < 0.1',
        total_df.filter(
            (pl.col('susie_cs') >= 0) & ~pl.col('finemap_pip').is_null()
        ).groupby([
            'phenotype', 'region', 'susie_cs'
        ]).agg(
            pl.col('finemap_pip').sum()
        ).with_column(
            (pl.col('finemap_pip') < .1).alias('discordant')
        ).groupby('discordant').agg(
            pl.count()
        ).with_column(
            pl.col('count').sum().alias('total')
        ).filter('discordant').select(
            pl.col('count')/pl.col('total')
        )
    )
    
    print(
        'Percentage of SuSiE CSes with FINEMAP PIP < 0.8',
        total_df.filter(
            (pl.col('susie_cs') >= 0) & ~pl.col('finemap_pip').is_null()
        ).groupby([
            'phenotype', 'region', 'susie_cs'
        ]).agg(
            pl.col('finemap_pip').sum()
        ).with_column(
            (pl.col('finemap_pip') < .8).alias('discordant')
        ).groupby('discordant').agg(
            pl.count()
        ).with_column(
            pl.col('count').sum().alias('total')
        ).filter('discordant').select(
            pl.col('count')/pl.col('total')
        )
    )

    # --- figures ----
    print('Plotting figures', flush=True)

    # Min abs corr across all CSes
    fig = bokeh.plotting.figure(
        width=1200,
        height=900,
        #title='SuSiE credible set purities',
        x_axis_label='Purity',
        y_axis_label='Number of credible sets',
        output_backend='svg'
    )
    fig.axis.axis_label_text_font_size = '36px'
    fig.title.text_font_size = '36px'
    fig.axis.major_label_text_font_size = '30px'
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color = None
    fig.toolbar_location = None
    step = 0.01
    left_edges = np.arange(0, 1 + step, step)
    ys = [np.sum((left_edge <= min_abs_corrs) & (min_abs_corrs < left_edge + step)) for left_edge in left_edges]
    ys[-2] += ys[-1]
    ys = ys[:-1]
    left_edges = left_edges[:-1]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+step)

    print('Exporting min abs corr plots', flush=True)
    bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/cs_min_abs_corrs.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/cs_min_abs_corrs.svg')

    fig_height=1200
    fig = bokeh.plotting.figure(
        width=fig_height,
        height=fig_height,
        #title='SuSiE alpha vs PIP',
        x_axis_label='Largest alpha',
        y_axis_label='PIP',
        output_backend='svg'
    )

    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color = None
    fig.toolbar_location = None
    fig.title.text_font_size = '36px'
    fig.axis.axis_label_text_font_size = '36px'
    fig.axis.major_label_text_font_size = '30px'
    bin_size = .02

    alpha_pip_df = total_df.filter(pl.col('susie_pip') >= 0.05)

    bins = bokeh.util.hex.hexbin(alpha_pip_df['susie_alpha'].to_numpy(), alpha_pip_df['susie_pip'].to_numpy(), size=bin_size)

    palette = [linear_int_interpolate((134,204,195), (9,41,46), i/254) for i in range(-1, 255)]
    cmap = bokeh.transform.log_cmap(
        'counts',
        palette = palette,
        low=1,
        high=max(bins.counts),
        low_color=(255, 255, 255)
    )
    color_mapper = bokeh.models.LogColorMapper(
        palette = palette,
        low=1,
        high=max(bins.counts)
    )

    fig.hex_tile(q='q', r='r', size=bin_size, line_color=None, source=bins, fill_color=cmap)
    color_bar = bokeh.models.ColorBar(
        title='Number of loci',
        color_mapper = color_mapper,
        width=35,
        height=fig_height//3,
        major_label_text_font_size = '30px',
        title_text_font_size='30px',
        ticker = bokeh.models.tickers.FixedTicker(ticks=[0, 10, 100, 1000, 10000, 100000]),
        formatter = bokeh.models.formatters.NumeralTickFormatter(format = '0a'),
        location='bottom_right'
    )
    fig.add_layout(color_bar, 'right')

    graphing_utils.resize(fig, 5000/1200, legend=False)
    bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/susie_alpha_v_pip.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/susie_alpha_v_pip.svg')

    fig = bokeh.plotting.figure(
        width=1200,
        height=900,
        #title='SuSiE total CP contributions',
        x_axis_label='CP',
        y_axis_label='Fraction of total CP contributions',
        output_backend='svg'
    )

    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color = None
    fig.toolbar_location = None
    fig.title.text_font_size = '36px'
    fig.axis.axis_label_text_font_size = '36px'
    fig.axis.major_label_text_font_size = '30px'

    step = 0.01
    left_edges = np.arange(0, 1 + step, step)
    total_susie_pip = total_df.filter(pl.col('susie_cs') >= 1)['susie_alpha'].sum()
    ys = [
        total_df.filter(
            (pl.col('p_val') <= 5e-8) &
            (pl.col('susie_cs') >= 1) &
            (float(left_edge) <= pl.col('susie_alpha')) &
            (pl.col('susie_alpha') < float(left_edge) + step)
        )['susie_alpha'].sum()/total_susie_pip
        for left_edge in left_edges
    ]
    ys[-2] += ys[-1]
    ys = ys[:-1]
    left_edges = left_edges[:-1]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+step)

    graphing_utils.resize(fig, 5000/1200, legend=False)
    bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/susie_alpha_histogram.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/susie_alpha_histogram.svg')

    fig = bokeh.plotting.figure(
        width=1200,
        height=900,
        #title='FINEMAP total CP contributions',
        x_axis_label='CP',
        y_axis_label='Fraction of total CP contributions',
        output_backend='svg'
    )

    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color = None
    fig.toolbar_location = None
    fig.title.text_font_size = '36px'
    fig.axis.axis_label_text_font_size = '36px'
    fig.axis.major_label_text_font_size = '30px'

    step = 0.01
    left_edges = np.arange(0, 1 + step, step)
    total_finemap_pip = total_df['finemap_pip'].sum()
    ys = [
        total_df.filter(
            (pl.col('p_val') <= 5e-8) &
            (float(left_edge) <= pl.col('finemap_pip')) &
            (pl.col('finemap_pip') < float(left_edge) + step)
        )['finemap_pip'].sum()/total_finemap_pip
        for left_edge in left_edges
    ]
    ys[-2] += ys[-1]
    ys = ys[:-1]
    left_edges = left_edges[:-1]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+step)

    graphing_utils.resize(fig, 5000/1200, legend=False)
    bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/finemap_pip_histogram.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/finemap_pip_histogram.svg')

    fig = bokeh.plotting.figure(
        width=1200,
        height=900,
        #title='FINEMAP total CPs for SuSiE credible sets',
        x_axis_label='FINEMAP total CP',
        y_axis_label='Credible set count',
        output_backend='svg'
    )
    fig.title.text_font_size = '36px'
    fig.axis.axis_label_text_font_size = '36px'
    fig.axis.major_label_text_font_size = '30px'
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.ygrid.grid_line_color = None
    fig.xgrid.grid_line_color = None
    fig.toolbar.logo = None
    fig.toolbar_location = None
    finemap_cs_coverages = total_df.filter(
        (pl.col('susie_cs') >= 0) & ~pl.col('finemap_pip').is_null()
    ).groupby([
        'phenotype', 'region', 'susie_cs'
    ]).agg(
        pl.col('finemap_pip').sum()
    )['finemap_pip'].to_numpy()
    print(f'Omitting {np.sum(finemap_cs_coverages >= 1 + step)} cses from this graph for having too high summed CSes')
    finemap_cs_coverages = finemap_cs_coverages[finemap_cs_coverages < 1 + step]
    max_coverage = np.max(finemap_cs_coverages)

    step = 0.01
    left_edges = np.arange(0, max_coverage + step, step)
    ys = [np.sum((left_edge <= finemap_cs_coverages) & (finemap_cs_coverages < left_edge + step)) for left_edge in left_edges]
    assert ys[-1] == 0
    ys = ys[:-1]
    left_edges = left_edges[:-1]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+step)

    print('Exporting FINEMAP CS PIP plots', flush=True)
    bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/susie_cs_finemap_total_pips.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/susie_cs_finemap_total_pips.svg')

    for var_type, filter_cond in [
        ('STR', pl.col('is_STR')),
        ('SNP/Indel', ~pl.col('is_STR'))
    ]:
        fig_height=1200
        fig = bokeh.plotting.figure(
            title=f'{var_type} CPs FINEMAP v SuSiE',
            width=fig_height,
            height=fig_height,
            y_axis_label = 'FINEMAP CP',
            x_axis_label = 'SuSiE CP',
            x_range=[0,1],
            y_range=[0,1],
        )
        fig.title.text_font_size = '36px'
        fig.axis.axis_label_text_font_size = '36px'
        
        min_text_size = '30px'
        fig.axis.major_label_text_font_size = min_text_size

        graph_df = total_df.with_columns([
            pl.when(
                pl.col('susie_cs') > 0
            ).then(
                pl.col('susie_alpha')
            ).otherwise(
                pl.lit(0)
            ).alias('susie_pip'),
            pl.max([pl.col('p_val'), 1e-300]).alias('p_val')
        ]).filter(
            (pl.col('p_val') <= 5e-8) &
            filter_cond & (
                (pl.col('finemap_pip') > 0.0) | (pl.col('susie_pip') > 0.0)
            )
        )

        n_rects = 100
        graph_df = graph_df.select([
            (pl.col('finemap_pip')*n_rects).floor()/n_rects,
            (pl.col('susie_pip')*n_rects).floor()/n_rects
        ]).groupby(['finemap_pip', 'susie_pip']).agg(pl.count())

        palette = [linear_int_interpolate((134,204,195), (9,41,46), i/254) for i in range(-1, 255)]
        cmap = bokeh.transform.log_cmap(
            'count',
            palette = palette,
            low=1,
            high=max(graph_df['count'].to_numpy()),
            low_color=(255, 255, 255)
        )
        color_mapper = bokeh.models.LogColorMapper(
            palette = palette,
            low=1,
            high=max(graph_df['count'].to_numpy())
        )

        cds = bokeh.models.ColumnDataSource(dict(
            left=graph_df['susie_pip'].to_numpy(), right=graph_df['susie_pip'].to_numpy() + 1/n_rects,
            bottom=graph_df['finemap_pip'].to_numpy(), top=graph_df['finemap_pip'].to_numpy() + 1/n_rects,
            count=graph_df['count'].to_numpy()
        ))
        fig.quad(
            left='left', right='right', bottom='bottom', top='top', source=cds, fill_color=cmap, line_width=0
        )

        thresh = .05
        both_df = graph_df.filter(
            (pl.col('finemap_pip') >= 1-thresh) &
            (pl.col('susie_pip') >= 1-thresh)
        )
        n_both = both_df['count'].sum()

        only_finemap_df = graph_df.filter(
            (pl.col('finemap_pip') >= 1-thresh) &
            (pl.col('susie_pip') <= thresh)
        )
        n_only_finemap = only_finemap_df['count'].sum()

        only_susie_df = graph_df.filter(
            (pl.col('susie_pip') >= 1-thresh) &
            (pl.col('finemap_pip') <= thresh)
        )
        n_only_susie = only_susie_df['count'].sum()

        line_width = 18
        xs = np.arange(0, 1, 0.0001)
        fig.line(
            xs,
            xs,
            line_width=line_width
        )

        fig.line(
            [1-thresh, 1-thresh],
            [1-thresh, 1],
            color='orange',
            width=line_width
        )
        fig.line(
            [1-thresh, 1],
            [1-thresh, 1-thresh],
            color='orange',
            width=line_width
        )

        fig.add_layout(bokeh.models.Title(
            text=f'# {var_type}s: {n_both}', align='right',
            text_font_size=min_text_size
        ), 'above')
        fig.line(
            [1-thresh, 1],
            [thresh, thresh],
            color='orange',
            width=line_width
        )
        fig.line(
            [1-thresh, 1-thresh],
            [0, thresh],
            color='orange',
            width=line_width
        )
        fig.add_layout(bokeh.models.Title(
            text=f'# {var_type}s: {n_only_susie}', align='right',
            text_font_size=min_text_size
        ), 'right')

        fig.line(
            [thresh, thresh],
            [1-thresh, 1],
            color='orange',
            width=line_width
        )
        fig.line(
            [0, thresh],
            [1-thresh, 1-thresh],
            color='orange',
            width=line_width
        )
        fig.add_layout(bokeh.models.Title(
            text=f'# {var_type}s: {n_only_finemap}', align='left',
            text_font_size=min_text_size
        ), 'above')

        color_bar = bokeh.models.ColorBar(
            title='Number of loci',
            color_mapper = color_mapper,
            width=35,
            height=fig_height//3,
            major_label_text_font_size = min_text_size,
            ticker = bokeh.models.tickers.FixedTicker(ticks=[0, 10, 100, 1000, 10000, 100000]),
            formatter = bokeh.models.formatters.NumeralTickFormatter(format = '0a'),
            title_text_font_size='30px',
            location='bottom_right'
        )
        fig.add_layout(color_bar, 'right')
        graphing_utils.resize(fig, 5000/1200, legend=False)
        fig.toolbar_location = None
        fig.background_fill_color = None
        fig.border_fill_color = None
        fig.grid.grid_line_color=None
        bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/consistency/finemap_v_susie_consistency_{var_type.split("/")[0]}.png')
        bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/consistency/finemap_v_susie_consistency_{var_type.split("/")[0]}.svg')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--putatively-causal-hits', action='store_true', default=False)
    parser.add_argument('--mpv', action='store_true', default=False)
    parser.add_argument('--first-pass', action='store_true', default=False)
    parser.add_argument('--regenerate', action='store_true', default=False)
    parser.add_argument('--phenotype')
    parser.add_argument('--no-assert', action='store_true', default=False)
    args = parser.parse_args()
    assert args.putatively_causal_hits + args.mpv + args.first_pass == 1

    if not args.regenerate:
        assert not args.no_assert

    assert args.regenerate == (args.phenotype is not None)

    if args.mpv:
        mpv_comparison()
    elif args.first_pass:
        if args.regenerate:
            first_pass_df(args.phenotype, not args.no_assert)
        else:
            first_pass_comparison()
    else:
        assert args.putatively_causal_hits
        if args.regenerate:
            putatively_causal_hits_df(args.phenotype, not args.no_assert)
        else:
            putatively_causal_hits_comparison()
