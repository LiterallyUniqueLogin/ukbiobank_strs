#!/usr/bin/env python3

import argparse
import glob
import os
import os.path
import re

import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

import phenotypes

ukb = os.environ['UKB']

corr_cutoff = .8
p_val_thresh = 5e-8

def load_susie(phenotype, results_regions_dir, colnames_regions_dir = None):
    dfs = []
    unconverged_regions = []
    underexplored_regions = []
    #unfinished_regions = []
    count = 0

    str_assocs = pl.scan_csv(
        f'{ukb}/association/results/{phenotype}/my_str/results.tab',
        sep='\t',
    ).select([
        pl.lit(phenotype).alias('phenotype'),
        'chrom',
        'pos',
        pl.col(f'p_{phenotype}').alias('p_val'),
        pl.lit(True).alias('is_STR'),
        pl.lit(None).cast(int).alias('reflen'),
        pl.lit(None).cast(int).alias('altlen')
    ])

    snp_assocs = pl.scan_csv(
        f'{ukb}/association/results/{phenotype}/plink_snp/results.tab',
        sep='\t',
        null_values='NA',
    ).select([
        pl.col('#CHROM').alias('chrom'),
        pl.col('POS').alias('pos'),
        pl.col('REF').str.lengths().cast(int).alias('reflen'),
        pl.col('ALT').str.lengths().cast(int).alias('altlen'),
        pl.col('P').alias('p_val'),
    ]).groupby(['chrom', 'pos', 'reflen', 'altlen']).agg([
        pl.col('p_val').min().alias('p_val'),
    ]).with_columns([
        pl.lit(phenotype).alias('phenotype'),
        pl.lit(False).alias('is_STR')
    ]).select(['phenotype', 'chrom', 'pos', 'p_val', 'is_STR', 'reflen', 'altlen'])

    assocs = pl.concat([str_assocs, snp_assocs])#.filter(pl.col('p_val') <= p_val_thresh)

    dirlist = os.listdir(results_regions_dir)
    for dir_ in dirlist:
        match = re.match('^([0-9]+)_([0-9]+)_([0-9]+)$', dir_)
        if not match:
            continue
        chrom, start, end = [int(match[idx]) for idx in (1,2,3)]
        if colnames_regions_dir:
            assert os.path.exists(f'{colnames_regions_dir}/{chrom}_{start}_{end}')
        if os.path.exists(f'{results_regions_dir}/{chrom}_{start}_{end}/no_strs'):
            continue

        converged_fname = f'{results_regions_dir}/{chrom}_{start}_{end}/converged.txt'
        if not os.path.exists(converged_fname):
            #unfinished_regions.append((phenotype, chrom, start, end))
            continue
        with open(converged_fname) as converged_file:
            if not next(converged_file).strip() == 'TRUE':
                unconverged_regions.append((phenotype, chrom, start, end))
                continue
        #print(f'Loading mean_platelet_volume region {chrom}:{start}-{end}', flush=True)
        if not colnames_regions_dir:
            colnames_fname = f'{results_regions_dir}/{chrom}_{start}_{end}/colnames.txt'
        else:
            colnames_fname = f'{colnames_regions_dir}/{chrom}_{start}_{end}/colnames.txt'
        '''
        if not os.path.exists(colnames_fname):
            colnames_fname = f'{colnames_fname}.normal_run'
        '''
        with open(colnames_fname) as var_file:
            susie_vars = [line.strip() for line in var_file if line.strip()]
        alphas = pl.scan_csv(
            f'{results_regions_dir}/{chrom}_{start}_{end}/alpha.tab',
            sep='\t',
            has_header=False
        ).collect().to_numpy().T
        n_alphas = alphas.shape[1]
        susie_pips=1-np.prod(1-alphas, axis=1)
        if not susie_pips.shape[0] == len(susie_vars):
            print(dir_, colnames_fname)
            assert False
        susie_idx = np.arange(len(susie_vars)) + 1
        susie_df = pl.DataFrame({
            'varname': susie_vars,
            'susie_pip': susie_pips,
            'susie_alpha': np.zeros(len(susie_vars)),
            'susie_cs': [-1]*len(susie_vars),
            'susie_idx': susie_idx,
            #**{ f'alpha_{i}': alphas[:, i] for i in range(n_alphas) }
        }).lazy()

        df = susie_df.with_columns([
            pl.col('varname').str.extract('^[^_]*_([^_]*)', 1).cast(int).alias('pos'),
            pl.col('varname').str.extract('^[^_]*_[^_]*_([^_]*)_.*', 1).str.lengths().cast(int).alias('reflen'),
            pl.col('varname').str.extract('^[^_]*_[^_]*_[^_]*_([^_]*)', 1).str.lengths().cast(int).alias('altlen'),
            pl.col('varname').str.contains('^STR').alias('is_STR'),
            pl.lit(f'{phenotype}_{chrom}_{start}_{end}').alias('region'),
            pl.lit(chrom).alias('chrom').cast(int),
            pl.lit(phenotype).alias('phenotype')
        ])#.sort('susie_idx')

        real_cs_count = 0
        for cs_fname in glob.glob(f'{results_regions_dir}/{chrom}_{start}_{end}/cs*.txt'):
            break
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
                underexplored_regions.append((phenotype, chrom, start, end))
            # could worry about variants being in multiple CSes
            df = df.with_column(
                pl.when(pl.col('susie_idx').is_in(cs_susie_idx))
                  .then(cs_id)
                  .otherwise(pl.col('susie_cs'))
                  .alias('susie_cs')
            )
        dfs.append(df)
        count += 1
        if count == 2:
            break
    dfs = [df.select(
        pl.col('*').exclude('^alpha.*$')
    ) for df in dfs]
    #print(pl.concat(dfs).filter((pl.col('varname') == 'SNP_25436941_T_C') | (pl.col('varname') == 'STR_25436639') | (pl.col('varname') == 'SNP_116396404_G_A') | (pl.col('varname') == 'STR_116545022')).collect().to_pandas())
    #print(assocs.filter(((pl.col('pos') == 25436941) & ~pl.col('is_STR')) | ((pl.col('pos') == 25436639) & pl.col('is_STR')) | ((pl.col('pos') == 116396404) & ~pl.col('is_STR')) | ((pl.col('pos') == 116545022) & pl.col('is_STR'))).collect().to_pandas())
    total_df = pl.concat(dfs).join(
        assocs,
        how='left',
        on=['phenotype', 'chrom', 'is_STR', 'pos', 'reflen', 'altlen']
    )#.collect()
    print(total_df.describe_optimized_plan())
    total_df = total_df.collect()
    print(total_df.to_pandas())#filter((pl.col('varname') == 'SNP_25436941_T_C') | (pl.col('varname') == 'STR_25436639') | (pl.col('varname') == 'SNP_116396404_G_A') | (pl.col('varname') == 'STR_116545022')).rechunk().to_pandas())
    exit()
    print('unconverged_regions: ', unconverged_regions)
    print('underexplored_regions: ', underexplored_regions)
    return total_df

def get_putatively_causal_regions():
    causal_df = pl.scan_csv(
        '{ukb}/post_finemappinng/results/validated/putatively_causal_hits.tab',
        sep='\t'
    )
    pheno_summaries = []
    for phenotype in phenotypes.phenotypes_in_use:
        pheno_summary = pl.scan_csv(
            f'{ukb}/finemapping/sumarry/{phenotype}_table.tab',
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
    ).select([
        'phenotype', 'signal_region'
    ]).drop_duplicates().collect().to_dict(False)

    phenos, regions = [
        causal_df[col] for col in ('phenotype', 'signal_region')
    ]
    return list(zip(phenos, regions))


'''
def putatively_causal_hits_comparison():
    putatively_causal_regions = get_putatively_causal_regions()
    for phenotype in phenotypes.phenotypes_in_use:
    .select([
        'phenotype',
    print("Loading baseline", flush=True)
    baseline = f'{ukb}/finemapping/susie_results/mean_platelet_volume_tol_0.001'
    original_susie = load_susie('mean_platelet_volume', original_susie_regions_dir, baseline).select([
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
        other_susie = load_susie('mean_platelet_volume', original_susie_regions_dir, other_susie_dir)
        assert 1 == other_susie.groupby(['phenotype', 'chrom', 'varname']).agg([pl.col('chrom').count().alias('count')]).sort('count')['count'].to_numpy()[-1]
        joined = other_susie.join(
            original_susie,
            how='outer',
            on = ['phenotype', 'chrom', 'varname']
        )
        assert 1 == joined.groupby(['phenotype', 'chrom', 'varname']).agg([pl.col('chrom').count().alias('count')]).sort('count')['count'].to_numpy()[-1]
        joined = joined.drop(['altlen', 'reflen', 'susie_idx', 'susie_pip', 'phenotype'])
        joined.to_csv(f'{ukb}/post_finemapping/intermediate_results/susie_{other_name}_comparison.tab', sep='\t')
'''

def mpv_comparison():
    original_susie_regions_dir = f'{ukb}/finemapping/susie_results/mean_platelet_volume'
    print("Loading baseline", flush=True)
    baseline = f'{ukb}/finemapping/susie_results/mean_platelet_volume_tol_0.001'
    original_susie = load_susie('mean_platelet_volume', baseline, original_susie_regions_dir).select([
        'phenotype', 'chrom', 'varname',
        pl.col('susie_alpha').alias('susie_baseline_alpha'),
        pl.col('susie_cs').alias('susie_baseline_cs')
    ])

    assert 1 == original_susie.groupby(['phenotype', 'chrom', 'varname']).agg([pl.col('chrom').count().alias('count')]).sort('count')['count'].to_numpy()[-1]

    for other_susie_dir, other_name in [
        (f'{ukb}/finemapping/susie_hardcall_results/mean_platelet_volume', 'hardcall'),
        #(f'{ukb}/finemapping/susie_results/mean_platelet_volume_tol_0.0001', 'tol_0.0001'),
        #(f'{ukb}/finemapping/susie_results/mean_platelet_volume_snp_str_ratio_1.5', 'ratio_1.5'),
        #(f'{ukb}/finemapping/susie_results/mean_platelet_volume_snp_str_ratio_4', 'ratio_4'),
        #(f'{ukb}/finemapping/susie_results/mean_platelet_volume_prior_var_0.2', 'prior_var_0.2'),
        #(f'{ukb}/finemapping/susie_results/mean_platelet_volume_prior_var_0.0005', 'prior_var_0.0005'),
        #(f'{ukb}/finemapping/susie_results/mean_platelet_volume_res_var_0.95', 'res_var_0.95'),
        #(f'{ukb}/finemapping/susie_results/mean_platelet_volume_res_var_0.8', 'res_var_0.8'),
    ]:
        print(f"Loading {other_name}", flush=True)
        other_susie = load_susie('mean_platelet_volume', other_susie_dir, original_susie_regions_dir)
        assert 1 == other_susie.groupby(['phenotype', 'chrom', 'varname']).agg([pl.col('chrom').count().alias('count')]).sort('count')['count'].to_numpy()[-1]
        joined = other_susie.join(
            original_susie,
            how='outer',
            on = ['phenotype', 'chrom', 'varname']
        )
        assert 1 == joined.groupby(['phenotype', 'chrom', 'varname']).agg([pl.col('chrom').count().alias('count')]).sort('count')['count'].to_numpy()[-1]
        joined = joined.drop(['altlen', 'reflen', 'susie_idx', 'susie_pip', 'phenotype'])
        joined.to_csv(f'{ukb}/post_finemapping/intermediate_results/susie_{other_name}_comparison.tab', sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--putatively-causal-hits', action='store_true', default=False)
    parser.add_argument('--mpv', action='store_true', default=False)
    args = parser.parse_args()
    assert args.putatively_causal_hits + args.mpv == 1
    if args.mpv:
        mpv_comparison()
    else:
        assert args.putatively_causal_hits
        putatively_causal_hits_comparison()
