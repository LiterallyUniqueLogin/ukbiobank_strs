#!/usr/bin/env python3

import argparse
import glob
import math
import os
import os.path
import pathlib
import time

import bokeh.models
import bokeh.io
import bokeh.layouts
import bokeh.plotting
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

corr_cutoff = .8
p_val_thresh = 5e-8

def all_regions(phenotype):
    return zip(*list(pl.read_csv(
        f'{ukb}/signals/regions/{phenotype}.tab',
        sep='\t'
    ).filter(
        pl.col('any_strs')
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

def load_susie(results_regions_dir, colnames_regions_dir = None, regions = None, phenotype=None, original=False, return_corrs = False):
    assert (phenotype is not None) + (regions is not None) == 1
    assert (colnames_regions_dir is not None) + original <= 1
    dfs = []
    unconverged_regions = []
    underexplored_regions = []
    unfinished_regions = []
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

        alphas = pl.scan_csv(
            f'{results_regions_dir}/{region}/alpha.tab',
            sep='\t',
            has_header=False
        ).collect().to_numpy().T
        n_alphas = alphas.shape[1]
        susie_pips=1-np.prod(1-alphas, axis=1)
        if not susie_pips.shape[0] == len(susie_vars):
            print(results_regions_dir, colnames_fname)
            assert False
        susie_idx = np.arange(len(susie_vars)) + 1
        susie_df = pl.DataFrame({
            'varname': susie_vars,
            #'susie_pip': susie_pips,
            'susie_alpha': np.zeros(len(susie_vars)),
            'susie_cs': [-1]*len(susie_vars),
            'susie_idx': susie_idx,
            **{ f'alpha_{i}': alphas[:, i] for i in range(n_alphas) }
        }).lazy()

        df = susie_df.with_columns([
            pl.lit(region).alias('region'),
            pl.lit(chrom).alias('chrom').cast(int),
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
    '''
    putatively_causal_regions = get_putatively_causal_regions()
    phenos_to_regions = {}

    for phenotype in phenotypes.phenotypes_in_use:
        phenos_to_regions[phenotype] = []

    for (phenotype, region, chrom) in putatively_causal_regions:
        phenos_to_regions[phenotype].append((region, chrom))

    pheno_dfs = []
    #for count, (phenotype, regions) in [(0, ('aspartate_aminotransferase', phenos_to_regions['aspartate_aminotransferase'])), (1, ('total_bilirubin', phenos_to_regions['total_bilirubin']))]:
    #for count, (phenotype, regions) in enumerate(phenos_to_regions.items()):
    #for count, (phenotype, regions) in [(0, ('mean_platelet_volume', phenos_to_regions['mean_platelet_volume']))]:
    #print(f"Loading phenotype #{count+1} ({phenotype})", flush=True)
        if len(regions) == 0:
            continue
    '''

    regions = get_putatively_causal_regions(phenotype)
    if len(regions) == 0:
        pathlib.Path(f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab').touch()
        pathlib.Path(f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab.empty').touch()
        return

    str_assocs = pl.scan_csv(
        f'{ukb}/association/results/{phenotype}/my_str/results.tab',
        sep='\t',
    ).select([
        'chrom',
        'pos',
        ('STR_' + pl.col('pos').cast(str)).alias('varname'),
        pl.lit(True).alias('is_STR'),
        pl.col(f'p_{phenotype}').alias('p_val'),
    ])

    other_ethnicity_assocs = None
    other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']
    for ethnicity in other_ethnicities:
        one_other_ethnicity = pl.scan_csv(
            f'{ukb}/association/results_finemapped_only/{ethnicity}/{phenotype}/my_str/results.tab',
            sep='\t',
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
    print('hardcall SuSiE')
    hardcall_susies = load_susie(f'{ukb}/finemapping/susie_hardcall_results/{phenotype}', regions=regions)
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
    print('gt threshold FINEMAP')
    gt_threshold_finemaps      = load_finemap(f'{ukb}/finemapping/finemap_results_threshold_0.0005/{phenotype}', regions=regions)

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
        gt_threshold_finemaps,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_gt_thresh'
    ).drop('region_gt_thresh').join(
        original_susies,
        how='outer',
        on=['chrom', 'varname'],
        suffix='_extra'
    ).drop('region_extra').join(
        hardcall_susies,
        how='outer',
        on = ['chrom', 'varname'],
        suffix='_hardcall'
    ).drop('region_hardcall').join(
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
    #pheno_dfs.append(pheno_df)

    #total_df = pl.concat(pheno_dfs).select([ # choose col order
    pheno_df.select([
        'phenotype',
        'chrom',
        'region',
        'pos',
        'is_STR',
        'varname',
        'p_val',
        'susie_cs',
        'susie_alpha',
        'susie_cs_hardcall',
        'susie_alpha_hardcall',
        'susie_cs_ratio',
        'susie_alpha_ratio',
        'finemap_pip',
        'finemap_pip_ratio',
        'finemap_pip_total_prob',
        'finemap_pip_prior_std_derived',
        'finemap_pip_prior_std_low',
        'finemap_pip_conv_tol',
        'finemap_pip_mac',
        'finemap_pip_gt_thresh',
        *[f'{ethnicity}_p_val' for ethnicity in other_ethnicities],
    ])
    pheno_df.to_csv(f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab', sep='\t')
    if should_assert:
        assert pheno_df.select((pl.col('region') == '').any().alias('region'))['region'].to_numpy()[0] == False
        assert np.all(~np.isnan(pheno_df['p_val'].to_numpy()))
        assert np.all(~np.isnan(pheno_df['finemap_pip'].to_numpy()))
        assert np.all(np.isnan(pheno_df['susie_alpha'].to_numpy()) == np.isnan(pheno_df['susie_alpha_hardcall'].to_numpy()))
        assert np.all(np.isnan(pheno_df['susie_alpha'].to_numpy()) == np.isnan(pheno_df['susie_alpha_ratio'].to_numpy()))
        assert np.all(1 == pheno_df.groupby(['phenotype', 'chrom', 'varname']).agg([pl.count()]).sort('count')['count'].to_numpy())
        for ethnicity in other_ethnicities:
            assert np.all(~np.isnan(pheno_df.filter('is_STR')[f'{ethnicity}_p_val'].to_numpy()))
    #return total_df

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
    ])
    print(' done', flush=True)

    total_df = total_df.filter(~pl.col('finemap_pip').is_null() & ~pl.col('susie_alpha').is_null())

    total_df = total_df.with_columns([
        pl.when(pl.col('susie_cs') > 0).then(pl.col('susie_alpha')).otherwise(0).alias('susie_alpha'),
        pl.when(pl.col('susie_cs_ratio') > 0).then(pl.col('susie_alpha_ratio')).otherwise(0).alias('susie_alpha_ratio'),
        pl.when(pl.col('susie_cs_hardcall') > 0).then(pl.col('susie_alpha_hardcall')).otherwise(0).alias('susie_alpha_hardcall'),
    ])
    total_df.filter(
        pl.col('is_STR') &
        (pl.col('p_val') <= 1e-10) &
        (pl.col('susie_alpha') >= .8) &
        (pl.col('finemap_pip') >= .8)
    ).drop('is_STR').to_csv(f'{ukb}/post_finemapping/intermediate_results/original_causal_STR_candidates.tab', sep='\t')

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
    ).drop(['susie_cs', 'susie_cs_hardcall', 'susie_cs_ratio', 'is_STR', 'varname'])
    pass_all_threshes.to_csv(f'{ukb}/post_finemapping/intermediate_results/concordant_causal_STR_candidates.tab', sep='\t')
    total_df.filter(
        pl.col('is_STR') &
        (pl.col('p_val') <= 1e-10) &
        (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols]) == 11)
    ).select(['phenotype', 'region', 'chrom', 'pos', 'p_val']).to_csv(f'{ukb}/post_finemapping/intermediate_results/strictly_concordant_causal_STR_candidates.tab', sep='\t')

    # plot replication numbers by finemapping category

    xs = []
    ys = []
    for name, condition in [
        ('either',
            pl.col('is_STR') &
            (pl.col('p_val') <= 1e-10) &
            ((pl.col('susie_alpha') >= .8) | (pl.col('finemap_pip') >= .8))
        ),
        ('both',
            pl.col('is_STR') &
            (pl.col('p_val') <= 1e-10) &
            ((pl.col('susie_alpha') >= .8) & (pl.col('finemap_pip') >= .8))
        ),
        ('candidate_STRs',
            pl.col('is_STR') &
            (pl.col('p_val') <= 1e-10) &
            (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols if 'ratio' not in col and 'prior_std_low' not in col]) == 8)
        ),
        ('strict candidates',
            pl.col('is_STR') &
            (pl.col('p_val') <= 1e-10) &
            (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols]) == 11)
        )
    ]:
        row = []
        for ethnicity in ['black', 'south_asian', 'chinese', 'irish', 'white_other']:
            xs.append((ethnicity, name))
            series = total_df.filter(condition).select(pl.col(f'{ethnicity}_p_val') <= 0.05).to_numpy()
            ys.append(np.sum(series)/series.shape[0])

        layout = bokeh.layouts.layout([
            bokeh.models.Title(title=f'correlations among {name} STRs'),
            row
        ])
        bokeh.io.export_png(layout, filename=f'{ukb}/post_finemapping/results/finemapping_{name}_ethnic_correlation.png')
    figure = bokeh.plotting.figure(
        width=1200,
        height=800,
        y_axis_label = 'Replication rate',
        x_range=bokeh.models.FactorRange(*xs),
        title='Replication rate by ethnicity and filtering',
        toolbar_location=None,
        tools=''
    )
    figure.vbar(x=xs, top=ys, width=0.9)
    figure.y_range.start=0
    figure.x_range.range_padding=0.1
    figure.xaxis.major_label_orientation = 1
    figure.grid.grid_line_color=None
    bokeh.io.export_png(figure, filename=f'{ukb}/post_finemapping/results/finemapping_results_ethnic_replication.png')

    '''
    # plot replication by p-value, broken down by finemapping category
    replication_data = total_df.select([
        'phenotype',
        'chrom',
        'pos',
        (-pl.max([pl.col('p_val'), 1e-300]).log10()).alias('-log10_p_val'),
        (pl.col('is_STR') & (pl.col('p_val') <= 1e-10) & ((pl.col('susie_alpha') >= .8) | (pl.col('finemap_pip') >= .8))).alias('singly_finemapped_STR'),
        (pl.col('is_STR') & (pl.col('p_val') <= 1e-10) & (pl.col('susie_alpha') >= .8) & (pl.col('finemap_pip') >= .8)).alias('doubly_finemapped_STR'),
        (pl.col('is_STR') & (pl.col('p_val') <= 1e-10) & (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols if 'ratio' not in col and 'prior_std_low' not in col]) == 8)).alias('candidate_STR'),
        (pl.col('is_STR') & (pl.col('p_val') <= 1e-10) & (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols]) == 11)).alias('strict_candidate_STR'),
        *[
            (pl.col(f'{ethnicity}_p_val') <= 0.05).alias(f'{ethnicity}_replication') for ethnicity in ['black', 'south_asian', 'chinese', 'irish', 'white_other']
        ],
    ]).filter('nearby_STR').to_csv(f'{ukb}/for_melissa.tab', sep='\t')
    exit()
    print(replication_data)
    print('plotting swarm plots ... ', end='', flush=True)
    f, axs = plt.subplots(1,5, figsize=(9,5), sharey=True)
    for count, col in enumerate(['nearby_STR', 'singly_finemapped_STR', 'doubly_finemapped_STR', 'candidate_STR', 'strict_candidate_STR']):
        print(f'{count} ', end='', flush=True)
        if count == 0: 
            ax = seaborn.violinplot(
                y='-log10_p_val',
                x='irish_replication',
                inner=None,
                ax=axs[count],
                data=replication_data.filter(col).to_pandas(),
                cut=0
            )
            continue
        #print(np.sum(replication_data[col].to_numpy()))
        
        ax = seaborn.swarmplot(
            y='-log10_p_val',
            x='irish_replication',
            dodge=True,
            ax=axs[count],
            data=replication_data.filter(col).to_pandas(),
            size=1
        )
        ax.set(title=col)
    plt.savefig(f'{ukb}/post_finemapping/results/finemapping_ethnic_replication_by_p_value_irish.png')
    print('done', flush=True)
    '''

    # plot fine-mapped SNPs, STRs broken down per pheno
    conditions = {
        'Variants fine-mapped by either the original SuSiE or FINEMAP runs':  ((
            (pl.col('p_val') <= 1e-10) &
            ((pl.col('susie_alpha') >= .8) | (pl.col('finemap_pip') >= .8))
        ), 'either'),
        'Variants fine-mapped by both the original SuSiE and FINEMAP runs': ((
            (pl.col('p_val') <= 1e-10) &
            (pl.col('susie_alpha') >= .8) &
            (pl.col('finemap_pip') >= .8)
        ), 'both'),
        'Variants fine-mapped under all selected fine-mapping conditions' : ((
            (pl.col('p_val') <= 1e-10) &
            (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols if 'ratio' not in col and 'prior_std_low' not in col]) == 8)
        ), 'concordant'),
        'Variants fine-mapped under all fine-mapping conditions': ((
            (pl.col('p_val') <= 1e-10) &
            (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols]) == 11)
        ), 'strictly_concordant')
    }

    '''
    pheno_pairs = [('biomarkers', pheno) if phenotypes.is_serum_biomarker(pheno) else ('blood cell counts', pheno) for pheno in phenotypes.phenotypes_in_use]
    sort_order_df = total_df.filter(conditions['Variants fine-mapped by either the original SuSiE or FINEMAP runs'][0]).groupby('phenotype').agg(pl.count())
    sort_key = lambda pheno_pair: (pheno_pair[0],
        sort_order_df.filter(pl.col('phenotype') == pheno_pair[1])['count'].to_numpy()[0] if sort_order_df.filter(pl.col('phenotype') == pheno_pair[1]).shape[0] > 0 else 0
    )
    sorted_pheno_pairs = sorted(pheno_pairs, key = sort_key)
    '''
    for title, (condition, file_title) in conditions.items():
        plot_df = total_df.filter(
            condition
        ).groupby([
            'phenotype', 'is_STR'
        ]).agg(pl.count())

        pheno_pairs = [('biomarkers', pheno) if phenotypes.is_serum_biomarker(pheno) else ('blood cell counts', pheno) for pheno in phenotypes.phenotypes_in_use]
        sort_order_df = total_df.filter(condition).groupby('phenotype').agg(pl.count())
        sort_key = lambda pheno_pair: (pheno_pair[0],
            sort_order_df.filter(pl.col('phenotype') == pheno_pair[1])['count'].to_numpy()[0] if sort_order_df.filter(pl.col('phenotype') == pheno_pair[1]).shape[0] > 0 else 0
        )
        sorted_pheno_pairs = sorted(pheno_pairs, key = sort_key)
        snps = np.array([
            plot_df.filter((pl.col('phenotype') == pheno) & ~pl.col('is_STR'))['count'].to_numpy()[0]
            if plot_df.filter((pl.col('phenotype') == pheno) & ~pl.col('is_STR')).shape[0] > 0 else 0
            for (_, pheno) in  sorted_pheno_pairs
        ])
        strs = np.array([
            plot_df.filter((pl.col('phenotype') == pheno) & pl.col('is_STR'))['count'].to_numpy()[0]
            if plot_df.filter((pl.col('phenotype') == pheno) & pl.col('is_STR')).shape[0] > 0 else 0
            for (_, pheno) in  sorted_pheno_pairs
        ])
        fig = bokeh.plotting.figure(
            x_range=bokeh.models.FactorRange(*sorted_pheno_pairs),
            title=title,
            tools='',
            x_axis_label = 'phenotype',
            y_axis_label = 'n fine-mapped variants',
            width=math.floor(4.25/2*1200),
            height=1200,
        )
        fig.toolbar_location = None
        fig.xaxis.major_label_orientation = 1.2
        fig.x_range.range_padding = 0.1

        fig.vbar(
            top=snps,
            x=sorted_pheno_pairs,
            legend_label='SNPs',
            color='#00B8FF',
            width=0.9
        )
        fig.vbar(
            top=strs+snps,
            bottom=snps,
            x=sorted_pheno_pairs,
            legend_label='STRs',
            color='#FF520D',
            width=0.9
        )
        fig.legend.location='top_left'
        fig.background_fill_color = None
        fig.border_fill_color = None
        fig.grid.grid_line_color = None
        fig.title.text_font_size = '30px'
        fig.axis.axis_label_text_font_size = '26px'
        fig.axis.major_label_text_font_size = '20px'
        fig.xaxis.group_text_font_size = '26px'
        fig.legend.label_text_font_size = '18px'

        graphing_utils.resize(fig, 5000/1200)
        bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/finemapping_results_by_pheno_{file_title}.png')

    #plot_upset_plots(total_df)
    
    # susie graphs
    for mapper, suffix, y_label in [
        ('SuSiE', 'ratio', '4x prior on SNPs'),
        ('SuSiE', 'hardcall', 'hardcall genotyping'),
        ('FINEMAP', 'ratio', '4x prior on SNPs'),
        ('FINEMAP', 'conv_tol', '10x stricter convergence tolerance'),
        ('FINEMAP', 'total_prob', 'assumption of 4 causal variants per region, not one'),
        ('FINEMAP', 'prior_std', '10x smaller assumption of effect size'),
        ('FINEMAP', 'mac', 'mac>=100 threshold'),
        ('FINEMAP', 'gt_thresh', 'p_val >= 5e-4 threshold instead of 5e-2')
    ]:
        if mapper == 'SuSiE':
            pip_col = 'susie_alpha'
            other_label = 'FINEMAP PIP'
            other_pip_col = 'finemap_pip'
        else:
            assert mapper == 'FINEMAP'
            pip_col = 'finemap_pip'
            other_label = 'SuSiE alpha'
            other_pip_col = 'susie_alpha'

        cs_STRs = total_df.filter(
            pl.col('is_STR') &
            (pl.col('p_val') <= 1e-10) &
            ((pl.col(pip_col) > 0) | (pl.col(f'{pip_col}_{suffix}') > 0))
        ).sort(pip_col).with_column(
            pl.max([pl.col('p_val'), 1e-300]).alias('p_val')
        )
        thresh = .025
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

        if n_new > 5:
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
        bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/consistency/{mapper.lower()}_consistency_{suffix}.png')
        bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/consistency/{mapper.lower()}_consistency_{suffix}.svg')

def first_pass_df(phenotype, should_assert):
    #pheno_dfs = []
    min_abs_corrs = []
    #for count, phenotype in [(0, 'aspartate_aminotransferase'), (1, 'total_bilirubin')]:
    #for count, phenotype in [(0, 'mean_platelet_volume')]:
    #for count, phenotype in enumerate(phenotypes.phenotypes_in_use):
    #    print(f"Loading phenotype #{count+1} ({phenotype})", flush=True)
    str_assocs = pl.scan_csv(
        f'{ukb}/association/results/{phenotype}/my_str/results.tab',
        sep='\t',
    ).select([
        'chrom',
        'pos',
        ('STR_' + pl.col('pos').cast(str)).alias('varname'),
        pl.lit(True).alias('is_STR'),
        pl.col(f'p_{phenotype}').alias('p_val'),
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
    ])
    assocs = pl.concat([snp_assocs, str_assocs])

    other_ethnicity_assocs = None
    other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']
    for ethnicity in other_ethnicities:
        one_other_ethnicity = pl.scan_csv(
            f'{ukb}/association/results_finemapped_only/{ethnicity}/{phenotype}/my_str/results.tab',
            sep='\t',
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
    #    pheno_dfs.append(pheno_df)

    #total_df = pl.concat(pheno_dfs).select([ # choose col order
    total_df = pheno_df.select([ # choose col order
        'phenotype',
        'chrom',
        'region',
        'pos',
        'is_STR',
        'varname',
        'p_val',
        'susie_cs',
        'susie_alpha',
        'finemap_pip',
        *[f'{ethnicity}_p_val' for ethnicity in other_ethnicities],
    ])
    total_df.to_csv(f'{ukb}/post_finemapping/intermediate_results/finemapping_all_concordance_{phenotype}.tab', sep='\t')
    #total_df.to_csv(f'{ukb}/post_finemapping/intermediate_results/finemapping_all_concordance.tab', sep='\t')
    np.save(f'{ukb}/post_finemapping/intermediate_results/susie_all_min_abs_corrs_{phenotype}.npy', np.array(min_abs_corrs))
    #np.save(f'{ukb}/post_finemapping/intermediate_results/susie_all_min_abs_corrs.npy', np.array(min_abs_corrs))
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
        ) for phenotype in phenotypes.phenotypes_in_use
    ])
    min_abs_corrs = np.concatenate([
        np.load(f'{ukb}/post_finemapping/intermediate_results/susie_all_min_abs_corrs_{phenotype}.npy')
        for phenotype in phenotypes.phenotypes_in_use
    ])

    '''
    # Min abs corr across all CSes
    fig = bokeh.plotting.figure(
        width=1200,
        height=1200,
        title='SuSiE credible set min absolute correlations',
        x_axis_label='min absolute correlation',
        y_axis_label='# credible sets',
    )
    fig.axis.axis_label_text_font_size = '30px'
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color = None
    fig.toolbar_location = None
    step = 0.01
    left_edges = np.arange(0, 1 + step, step)
    ys = [np.sum((left_edge <= min_abs_corrs) & (min_abs_corrs < left_edge + step)) for left_edge in left_edges]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+step)

    print('Exporting min abs corr plots', flush=True)
    bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/cs_min_abs_corrs.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/cs_min_abs_corrs.svg')

    # STR fraction per CS
    fraction_strs_df = total_df.filter(
        pl.col('susie_cs') >= 0
    ).groupby([
        'phenotype',
        'region',
        'susie_cs',
        'is_STR'
    ]).agg(
        pl.col('susie_alpha').sum()
    ).with_column(
        pl.col('susie_alpha').sum().over(['phenotype', 'region', 'susie_cs']).alias('total_alpha')
    )
    assert np.all(fraction_strs_df['total_alpha'].to_numpy() >= .9)
    joined_fraction_df = fraction_strs_df.filter(~pl.col('is_STR')).join(
        fraction_strs_df.filter(pl.col('is_STR')),
        how='outer',
        on=['phenotype', 'region', 'susie_cs'],
        suffix='_str'
    )
        
    assert np.all(joined_fraction_df.filter(
        ~pl.col('total_alpha').is_null() & ~pl.col('total_alpha_str').is_null()
    ).select(
        (pl.col('total_alpha') == pl.col('total_alpha_str')).alias('out')
    )['out'].to_numpy())

    fraction_strs = joined_fraction_df.select(
        pl.when(
            ~pl.col('total_alpha_str').is_null()
        ).then(
            pl.col('susie_alpha_str')/pl.col('total_alpha_str')
        ).otherwise(
            0
        ).alias('out')
    )['out'].to_numpy()

    pdf = scipy.stats.kde.gaussian_kde(fraction_strs)

    fig = bokeh.plotting.figure(
        height=1200,
        width=900,
        x_range = (0,1),
        toolbar_location=None
    )

    xs = np.arange(0, 1, 0.001)
    ys = pdf(xs)
    fig.varea(xs, np.zeros(xs.shape), ys, fill_color='orange')
    for count, perc in enumerate((.8, .9, .95)):
        q = np.quantile(fraction_strs, perc)
        fig.line([q, q], [0, max(ys)], color='black')
        fig.add_layout(bokeh.models.Label(
            x=q, y=max(ys)-3-count, text=f'{perc*100}% quantile'
        ))
    perc = np.sum(fraction_strs >= .95)/fraction_strs.shape[0]
    fig.line([.95, .95], [0, max(ys)], color='black')
    fig.add_layout(bokeh.models.Label(
        x=.95, y=max(ys)-3, text=f'{(1-perc)*100:.4f}% quantile'
    ))

    bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/STR_fraction_per_cs.png')
    '''
    exit()

    fig = bokeh.plotting.figure(
        width=1200,
        height=1200,
        title=f'FINEMAP total PIPs for SuSiE CSes with min_abs_corr >= {corr_cutoff}',
        x_axis_label='FINEMAP PIPs',
        y_axis_label='# credible sets',
    )
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.ygrid.grid_line_color = None
    fig.xgrid.grid_line_color = None
    fig.toolbar.logo = None
    fig.toolbar_location = None
    finemap_cs_coverages = total_df.filter(
        pl.col('susie_cs') >= 0
    ).groupby([
        'phenotype', 'region', 'susie_cs'
    ]).agg(
        pl.col('finemap_pip').sum()
    )['finemap_pip'].to_numpy()
    max_coverage = np.max(finemap_cs_coverages)

    step = 0.01
    left_edges = np.arange(0, max_coverage + step, step)
    ys = [np.sum((left_edge <= finemap_cs_coverages) & (finemap_cs_coverages < left_edge + step)) for left_edge in left_edges]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+step)

    print('Exporting FINEMAP CS PIP plots', flush=True)
    bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/susie_cs_finemap_total_pips.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/susie_cs_finemap_total_pips.svg')

    '''
    # susie pip vs alpha
    figure, ax = plt.subplots(1, 1, figsize=(12,12))

    if vartype == '':
        title_inset = 'all'
    else:
        title_inset = vartype
    if not all_:
        title_inset += f' credible-set-var (min_abs_corr >= {corr_cutoff})'
    ax.set_title(f'SuSiE alpha vs PIP {title_inset}')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('SuSiE alpha')
    ax.set_ylabel('SuSiE PIPs')
    idx = total_df['varname'].str.contains(f'^{vartype}')
    idx = idx & (total_df['susie_cs'] | all_)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(name='meh', colors=['#e5f5f9', '#99d8c9', '#2ca25f'])
    hexbin = ax.hexbin(total_df[idx, 'susie_alpha'].to_numpy().flatten(), total_df[idx, 'susie_pip'].to_numpy().flatten(), gridsize=100, cmap=cmap, bins='log')
    cb = figure.colorbar(hexbin, ax=ax)
    cb.set_label('counts')

    if vartype != '':
        fname_suffix = '_' + vartype
    else:
        fname_suffix = ''
    if not all_:
        fname_suffix += '_credible_set_vars'
    out_fname = f'{ukb}/export_scripts/results/susie_alpha_v_pip{fname_suffix}'
    print(f'Exporting var plot {out_fname}.svg', flush=True)
    plt.savefig(f'{out_fname}.svg')
    print(f'Exporting var plot {out_fname}.png', flush=True)
    plt.savefig(f'{out_fname}.png')
    '''

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
