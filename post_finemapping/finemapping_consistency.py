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

def load_susie(phenotype, results_regions_dir, colnames_regions_dir = None, regions = None):
    dfs = []
    unconverged_regions = []
    underexplored_regions = []
    unfinished_regions = []

    str_assocs = pl.scan_csv(
        f'{ukb}/association/results/{phenotype}/my_str/results.tab',
        sep='\t',
    ).select([
        pl.lit(phenotype).alias('phenotype'),
        'chrom',
        'pos',
        pl.col(f'p_{phenotype}').alias('p_val'),
        ('STR_' + pl.col('pos').cast(str)).alias('varname'),
        pl.lit(True).alias('is_STR')
    ])

    snp_assocs = pl.scan_csv(
        f'{ukb}/association/results/{phenotype}/plink_snp/results.tab',
        sep='\t',
        null_values='NA',
    ).select([
        pl.col('#CHROM').alias('chrom'),
        pl.col('POS').alias('pos'),
        ('SNP_' + pl.col('POS').cast(str) + '_' + pl.col('REF') + '_' +  pl.col('ALT')).alias('varname'),
        pl.col('P').alias('p_val'),
    ]).with_columns([
        pl.lit(phenotype).alias('phenotype'),
        pl.lit(False).alias('is_STR')
    ]).select(['phenotype', 'chrom', 'pos', 'p_val', 'varname', 'is_STR'])

    assocs = pl.concat([str_assocs, snp_assocs])

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
            unfinished_regions.append((phenotype, region))
            continue
        with open(converged_fname) as converged_file:
            if not next(converged_file).strip() == 'TRUE':
                unconverged_regions.append((phenotype, region))
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
            print(dir_, colnames_fname)
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
            pl.col('varname').str.extract('^[^_]*_([^_]*)', 1).cast(int).alias('pos'),
            pl.lit(f'{phenotype}_{region}').alias('region'),
            pl.lit(chrom).alias('chrom').cast(int),
            pl.lit(phenotype).alias('phenotype')
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
                underexplored_regions.append((phenotype, region))
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
    total_df = pl.concat(dfs).join(
        assocs,
        how='left',
        on=['phenotype', 'chrom', 'varname']
    ).drop('pos_right').collect()
    assert np.all(~np.isnan(total_df['p_val'].to_numpy()))
    return total_df

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
    ).select([
        'phenotype', 'signal_region', 'chrom'
    ]).drop_duplicates().collect().to_dict(False)

    phenos, regions, chroms = [
        causal_df[col] for col in ('phenotype', 'signal_region', 'chrom')
    ]
    return list(zip(phenos, regions, chroms))


def putatively_causal_hits_comparison():
    putatively_causal_regions = get_putatively_causal_regions()
    phenos_to_regions = {}
    for pheno in phenotypes.phenotypes_in_use:
        phenos_to_regions[pheno] = []
    for (pheno, region, chrom) in putatively_causal_regions:
        phenos_to_regions[pheno].append((region, chrom))
    original_susies = []
    hardcall_susies = []
    #for count, (pheno, regions) in [(0, ('alanine_aminotransferase', phenos_to_regions['alanine_aminotransferase'])), (1, ('total_bilirubin', phenos_to_regions['total_bilirubin']))]:
    for count, (pheno, regions) in enumerate(phenos_to_regions.items()):
        if len(regions) == 0:
            continue
        print(f"Loading phenotype #{count} ({pheno})", flush=True)
        original_susies.append(load_susie(pheno, f'{ukb}/finemapping/susie_results/{pheno}', regions=regions))
        hardcall_susies.append(load_susie(pheno, f'{ukb}/finemapping/susie_hardcall_results/{pheno}', regions=regions))

    original_susie = pl.concat(original_susies)
    hardcall_susie = pl.concat(hardcall_susies)
        
    assert 1 == original_susie.groupby(['phenotype', 'chrom', 'varname']).agg([pl.col('chrom').count().alias('count')]).sort('count')['count'].to_numpy()[-1]
    assert 1 == hardcall_susie.groupby(['phenotype', 'chrom', 'varname']).agg([pl.col('chrom').count().alias('count')]).sort('count')['count'].to_numpy()[-1]

    joined = hardcall_susie.join(
        original_susie,
        how='outer',
        on = ['phenotype', 'chrom', 'varname']
    )
    assert 1 == joined.groupby(['phenotype', 'chrom', 'varname']).agg([pl.col('chrom').count().alias('count')]).sort('count')['count'].to_numpy()[-1]
    joined = joined.drop(['susie_idx', 'susie_pip'])
    joined.to_csv(f'{ukb}/post_finemapping/intermediate_results/susie_putative_hardcall_comparison.tab', sep='\t')

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
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_tol_0.0001', 'tol_0.0001'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_snp_str_ratio_1.5', 'ratio_1.5'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_snp_str_ratio_4', 'ratio_4'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_prior_var_0.2', 'prior_var_0.2'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_prior_var_0.0005', 'prior_var_0.0005'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_res_var_0.95', 'res_var_0.95'),
        (f'{ukb}/finemapping/susie_results/mean_platelet_volume_res_var_0.8', 'res_var_0.8'),
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
        joined = joined.drop(['susie_idx', 'susie_pip', 'phenotype'])
        joined.to_csv(f'{ukb}/post_finemapping/intermediate_results/susie_mpv_{other_name}_comparison.tab', sep='\t')

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
