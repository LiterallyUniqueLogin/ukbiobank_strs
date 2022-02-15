#!/usr/bin/env python3

import argparse
import glob
import os
import os.path
import re
import time

import bokeh.models
import bokeh.io
import bokeh.plotting
import numpy as np
import polars as pl
import scipy.interpolate

import phenotypes

ukb = os.environ['UKB']

corr_cutoff = .8
p_val_thresh = 5e-8

def load_susie(results_regions_dir, colnames_regions_dir = None, regions = None):
    dfs = []
    unconverged_regions = []
    underexplored_regions = []
    unfinished_regions = []

    if regions is None:
        regions = []
        dirlist = os.listdir(results_regions_dir)
        for dir_ in dirlist:
            match = re.match('^([0-9]+)_[0-9]+_[0-9]+$', dir_)
            if not match:
                continue
            regions.append((dir_, match[1]))
    for (region, chrom) in regions:
        if colnames_regions_dir:
            assert os.path.exists(f'{colnames_regions_dir}/{region}')
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
        #print(f'Loading mean_platelet_volume region {chrom}:{start}-{end}', flush=True)
        if not colnames_regions_dir:
            colnames_fname = f'{results_regions_dir}/{region}/colnames.txt'
        else:
            colnames_fname = f'{colnames_regions_dir}/{region}/colnames.txt'
        if not os.path.exists(colnames_fname):
            colnames_fname = f'{colnames_fname}.normal_run'
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
    return pl.concat(dfs).drop(['susie_idx'])

def get_putatively_causal_regions():
    causal_df = pl.scan_csv(
        f'{ukb}/post_finemapping/results/validated/putatively_causal_STRs.tab',
        sep='\t'
    )
    pheno_summaries = []
    for phenotype in phenotypes.phenotypes_in_use:
        pheno_summary = pl.scan_csv(
            f'{ukb}/finemapping/summary/{phenotype}_table.tab',
            sep='\t'
        ).with_column(
            pl.lit(phenotype).alias('phenotype')
        ).select([
            'phenotype', 'chrom', 'start_pos', 'signal_region'
        ])
        pheno_summaries.append(pheno_summary)
    causal_df = causal_df.join(
        pl.concat(pheno_summaries),
        how='left',
        on=['phenotype', 'chrom', 'start_pos']
    ).filter(~(
        (
            (pl.col('phenotype') == 'total_bilirubin') &
            (pl.col('signal_region') == '12_19976272_22524428')
        ) |
        (
            (pl.col('phenotype') == 'mean_platelet_volume') &
            (pl.col('signal_region') == '17_2341352_2710113')
        )
    )).select([
        'phenotype', 'signal_region', 'chrom'
    ]).distinct().collect().to_dict(False)

    phenos, regions, chroms = [
        causal_df[col] for col in ('phenotype', 'signal_region', 'chrom')
    ]
    return list(zip(phenos, regions, chroms))

def putatively_causal_hits_df():
    putatively_causal_regions = get_putatively_causal_regions()
    phenos_to_regions = {}

    for phenotype in phenotypes.phenotypes_in_use:
        phenos_to_regions[phenotype] = []

    for (phenotype, region, chrom) in putatively_causal_regions:
        phenos_to_regions[phenotype].append((region, chrom))

    pheno_dfs = []
    #for count, (phenotype, regions) in [(0, ('mean_platelet_volume', phenos_to_regions['mean_platelet_volume']))]:
    #for count, (phenotype, regions) in [(0, ('alanine_aminotransferase', phenos_to_regions['alanine_aminotransferase'])), (1, ('total_bilirubin', phenos_to_regions['total_bilirubin']))]:
    for count, (phenotype, regions) in enumerate(phenos_to_regions.items()):
        print(f"Loading phenotype #{count+1} ({phenotype})", flush=True)
        if len(regions) == 0:
            continue
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

        original_susies = load_susie(f'{ukb}/finemapping/susie_results/{phenotype}', regions=regions)
        hardcall_susies = load_susie(f'{ukb}/finemapping/susie_hardcall_results/{phenotype}', regions=regions)
        ratio_susies = load_susie(
            f'{ukb}/finemapping/susie_results/{phenotype}_snp_str_ratio_4',
            colnames_regions_dir=f'{ukb}/finemapping/susie_results/{phenotype}',
            regions=regions
        )

        print('Collecting ... ', end='', flush=True)
        start = time.time()
        pheno_df = original_susies.join(
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
        ).with_column(
            pl.lit(phenotype).alias('phenotype')
        ).collect()
        print(f'done. Time: {(time.time() - start)/60:.2f}m', flush=True)
        pheno_dfs.append(pheno_df)

    total_df = pl.concat(pheno_dfs).select([ # choose col order
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
    ])
    assert np.all(~np.isnan(total_df['p_val'].to_numpy()))
    assert np.all(~np.isnan(total_df['susie_alpha'].to_numpy()))
    assert np.all(~np.isnan(total_df['susie_alpha_hardcall'].to_numpy()))
    assert np.all(1 == total_df.groupby(['phenotype', 'chrom', 'varname']).agg([pl.count()]).sort('count')['count'].to_numpy())
    total_df.to_csv(f'{ukb}/post_finemapping/intermediate_results/susie_putative_hardcall_comparison.tab', sep='\t')
    return total_df

def putatively_causal_hits_comparison(regenerate):
    if regenerate:
        total_df = putatively_causal_hits_df()
    else:
        total_df = pl.read_csv(f'{ukb}/post_finemapping/intermediate_results/susie_putative_hardcall_comparison.tab', sep='\t')

    pc_STRs = total_df.filter(
        pl.col('is_STR') &
        (pl.col('p_val') <= 1e-10) &
        (pl.col('susie_cs') >= 0) &
        (pl.col('susie_alpha') >= 0.8)
    )
    print(f'SuSiE putatively causal STRs: {pc_STRs.shape[0]}')
    hard_rep = pc_STRs.filter(
        (pl.col('susie_cs_hardcall') >= 0) &
        (pl.col('susie_alpha_hardcall') >= .8)
    )
    print(f'That replicate through hardcalls: {hard_rep.shape[0]}')
    ratio_rep = pc_STRs.filter(
        (pl.col('susie_cs_ratio') >= 0) &
        (pl.col('susie_alpha_ratio') >= .8)
    )
    print(f'That replicate through ratios: {ratio_rep.shape[0]}')

    
    for suffix, y_label in ('ratio', '4x prior on SNPs'), ('hardcall', 'hardcall genotyping'):
        cs_STRs = total_df.filter(
            pl.col('is_STR') &
            (pl.col('p_val') <= 1e-10) &
            ((pl.col('susie_cs') >= 0) | (pl.col(f'susie_cs_{suffix}') >= 0))
        ).sort('susie_alpha').with_column(
            pl.max([pl.col('p_val'), 1e-300]).alias('p_val')
        )
        thresh = .025
        n_both = cs_STRs.filter(
            (pl.col('susie_alpha') >= 1-thresh) &
            (pl.col(f'susie_alpha_{suffix}') >= 1-thresh)
        ).shape[0]

        n_not_rep = cs_STRs.filter(
            (pl.col('susie_alpha') >= 1-thresh) &
            (pl.col(f'susie_alpha_{suffix}') <= thresh)
        ).shape[0]

        n_new = cs_STRs.filter(
            (pl.col('susie_alpha') <= thresh) &
            (pl.col(f'susie_alpha_{suffix}') >= 1-thresh)
        ).shape[0]

        fig = bokeh.plotting.figure(
            width=1200,
            height=1200,
            y_axis_label = f'PIP under {y_label}',
            x_axis_label = 'original PIP',
            x_range=[0,1],
            y_range=[0,1]
        )
        fig.circle(
            cs_STRs['susie_alpha'].to_numpy(),
            cs_STRs[f'susie_alpha_{suffix}'].to_numpy(),
            size = -np.log10(cs_STRs['p_val'].to_numpy())/7.5,
            alpha = 0.25
        )
        spline = scipy.interpolate.UnivariateSpline(
            cs_STRs['susie_alpha'].to_numpy(),
            cs_STRs[f'susie_alpha_{suffix}'].to_numpy()
        )
        xs = np.arange(0, 1, 0.0001)
        fig.line(
            xs,
            spline(xs),
            #legend_label='approximation'
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
            text=f'# STRs: {n_both}', align='right'
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
            text=f'# STRs: {n_not_rep}', align='right'
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
                text=f'# STRs: {n_new}', align='left'
            ), 'above')
            
        fig.toolbar_location = None
        fig.background_fill_color = None
        fig.border_fill_color = None
        bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/susie_consistency_{suffix}.png')

def mpv_comparison():
    original_susie_regions_dir = f'{ukb}/finemapping/susie_results/mean_platelet_volume'
    print("Loading baseline", flush=True)
    baseline = f'{ukb}/finemapping/susie_results/mean_platelet_volume_tol_0.001'
    original_susie = load_susie(baseline, original_susie_regions_dir).select([
        'phenotype', 'chrom', 'varname',
        pl.col('susie_alpha').alias('susie_baseline_alpha'),
        pl.col('susie_cs').alias('susie_baseline_cs')
    ])

    assert 1 == original_susie.groupby(['phenotype', 'chrom', 'varname']).agg([pl.col('chrom').count().alias('count')]).sort('count')['count'].to_numpy()[-1]

    for other_susie_dir, other_name in [
        (f'{ukb}/finemapping/susie_hardcall_results/mean_platelet_volume', 'hardcall'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_tol_0.0001', 'tol_0.0001'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_snp_str_ratio_1.5', 'ratio_1.5'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_snp_str_ratio_4', 'ratio_4'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_prior_var_0.2', 'prior_var_0.2'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_prior_var_0.0005', 'prior_var_0.0005'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_res_var_0.95', 'res_var_0.95'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_res_var_0.8', 'res_var_0.8'),
    ]:
        print(f"Loading {other_name}", flush=True)
        other_susie = load_susie(other_susie_dir, original_susie_regions_dir)
        assert 1 == other_susie.groupby(['phenotype', 'chrom', 'varname']).agg([pl.col('chrom').count().alias('count')]).sort('count')['count'].to_numpy()[-1]
        joined = other_susie.join(
            original_susie,
            how='outer',
            on = ['phenotype', 'chrom', 'varname']
        )
        assert 1 == joined.groupby(['phenotype', 'chrom', 'varname']).agg([pl.col('chrom').count().alias('count')]).sort('count')['count'].to_numpy()[-1]
        joined = joined.drop(['susie_idx', 'susie_pip', 'phenotype'])
        joined.to_csv(f'{ukb}/post_finemapping/intermediate_results/susie_mpv_{other_name}_comparison.tab', sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--putatively-causal-hits', action='store_true', default=False)
    parser.add_argument('--mpv', action='store_true', default=False)
    parser.add_argument('--regenerate', action='store_true', default=False)
    args = parser.parse_args()
    assert args.putatively_causal_hits + args.mpv == 1
    if args.mpv:
        mpv_comparison()
    else:
        assert args.putatively_causal_hits
        putatively_causal_hits_comparison(regenerate=args.regenerate)
