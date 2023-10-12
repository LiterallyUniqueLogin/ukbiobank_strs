#!/usr/bin/env python3

import argparse

import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('simulations_tsv')
args = parser.parse_args()

def load_finemap(snp_file, log_sss, creds):
    with open(log_sss) as log:
        found_n_causal = False
        for line in log:
            if 'n-causal' not in line:
                continue
            found_n_causal = True
            n_causal = int(line.split()[-1])
            with open(snp_file, 'rb') as f:
                n_vars = sum(1 for _ in f) - 1
            if n_causal == n_vars:
                # not underexplored if there are no more variants to explore
                break
            if any(cred.endswith(f'.cred{n_causal}') for cred in creds):
                print(f'Underexplored region! {snp_file}', flush=True)
            break
        assert found_n_causal

    return pl.scan_csv(
        snp_file,
        sep=' '
    ).select([
        pl.col('rsid').alias('varname'),
        pl.col('prob').cast(float).alias('finemap_pip'),
    ])

def load_susie(converged, colnames, alpha, V, CSes):
    if converged[-4:] == 'null':
        print('Missing region!', flush=True)
        return None
    with open(converged) as converged_file:
        if not next(converged_file).strip() == 'TRUE':
            print(f'Unconverged region! {converged}', flush=True)
            return None

    with open(colnames) as var_file:
        susie_vars = [line.strip() for line in var_file if line.strip()]

    alphas = np.genfromtxt(
        alpha,
        delimiter='\t',
    ).T
    if len(alphas.shape) == 0:
        alphas = alphas.reshape(-1, 1) # only one credible set
    estimated_signal_vars = np.genfromtxt(
        V,
        delimiter='\t'
    )
    if len(estimated_signal_vars.shape) == 0:
        estimated_signal_vars = np.array([estimated_signal_vars])

    n_alphas = alphas.shape[1]
    susie_idx = np.arange(len(susie_vars)) + 1
    df = pl.DataFrame({
        'varname': susie_vars,
        'susie_alpha': np.zeros(len(susie_vars)),
        'susie_cs': [-1]*len(susie_vars),
        'susie_idx': susie_idx,
        **{ f'alpha_{i}': alphas[:, i] for i in range(n_alphas) }
    }).lazy().sort('susie_idx')

    real_cs_count = 0
    for cs_fname in CSes:
        cs_id = int(cs_fname.split('cs')[-1].split('.')[0])
        with open(cs_fname) as cs_file:
            # susie uses 1 based indexing, python uses 0
            # make sure cs idxs are in increasing order
            cs_susie_idx = np.array([int(idx) for idx in next(cs_file).strip().split()])
            assert np.all(cs_susie_idx[1:] - cs_susie_idx[:-1] > 0)
            cs_susie_idx = pl.Series('cs_susie_idx', cs_susie_idx)
            next(cs_file) # skip cs credibility
            min_abs_corr, _, _ = [float(idx) for idx in next(cs_file).strip().split()]
        if min_abs_corr < corr_cutoff:
            continue
        df = df.with_column(
            pl.when(pl.col('susie_idx').is_in(cs_susie_idx))
              .then(pl.when(pl.col(f'alpha_{cs_id-1}') > pl.col('susie_alpha'))
                      .then(pl.col(f'alpha_{cs_id-1}'))
                      .otherwise(pl.col('susie_alpha')))
              .otherwise(pl.col('susie_alpha'))
              .alias('susie_alpha')
        )
        if estimated_signal_vars[cs_id-1] <= 1e-9:
            assert False, f'CS {cs_id} in region {colnames} has a pure CS with negligible signal, exiting'
        real_cs_count += 1
        if real_cs_count == 50:
            print(f'Underexplored region! {colnames}', flush=True)
        # could worry about variants being in multiple CSes
        df = df.with_column(
            pl.when(pl.col('susie_idx').is_in(cs_susie_idx))
              .then(cs_id)
              .otherwise(pl.col('susie_cs'))
              .alias('susie_cs')
        )
    return df.select(
        pl.col('*').exclude('^alpha.*$')
    ).drop(['susie_idx'])

corr_cutoff = 0.8

df = pl.read_csv(args.simulations_tsv, sep='\t')
finemap_dfs = []
susie_dfs = []
snp_assoc_dfs = []
str_assoc_dfs = []

# temporary workaround
def fix_file(unfixed_file):
    return (line for line in unfixed_file if line.strip())

for row in range(df.shape[0]):
    print(f'Loading simulation {row+1}/{df.shape[0]}', flush=True)

    with open(df['causal_vars_and_betas'][row]) as causal_file_unfixed:
        causal_file = fix_file(causal_file_unfixed)
        vars_ = next(causal_file).strip()
        betas = next(causal_file).strip()
        vars_betas = ' '.join([f'{var}:{beta}' for var, beta in zip(vars_.split(), betas.split())])
        done = False
        try:
            next(causal_file)
        except StopIteration:
            done = True
        assert done

    finemap_dfs.append(load_finemap(
        snp_file = df['snp_file'][row],
        log_sss = df['log_sss'][row],
        creds = df['creds'][row].split(','),
    ).with_columns([
        pl.lit(df['method'][row]).alias('method'),
        pl.lit(df['chrom'][row]).alias('chrom').cast(int),
        pl.lit(df['region'][row]).alias('region'),
        pl.lit(df['replicate'][row]).alias('replicate').cast(int),
        pl.lit(vars_).alias('causal_vars'),
        pl.lit(betas).alias('causal_betas'),
        pl.lit(vars_betas).alias('causal_vars_and_betas'),
    ]))

    one_susie_df = load_susie(
        df['converged'][row],
        colnames = df['colnames'][row],
        alpha = df['alpha'][row],
        V = df['V'][row],
        CSes = df['CSes'][row].split(',') if df['CSes'][row] is not None else [],
    )
    if one_susie_df is not None:
        one_susie_df = one_susie_df.with_columns([
            pl.lit(df['method'][row]).alias('method'),
            pl.lit(df['chrom'][row]).alias('chrom').cast(int),
            pl.lit(df['region'][row]).alias('region'),
            pl.lit(df['replicate'][row]).alias('replicate').cast(int),
            pl.lit(vars_).alias('causal_vars'),
            pl.lit(betas).alias('causal_betas'),
            pl.lit(vars_betas).alias('causal_vars_and_betas'),
        ])
        susie_dfs.append(one_susie_df)
    else:
        susie_dfs.append(susie_dfs[-1].collect()[:0, :].lazy())

    snp_assoc_dfs.append(pl.scan_csv(
        df['snp_assoc'][row],
        sep='\t',
        null_values = 'NA'
    ).select([
        pl.col('#CHROM').alias('chrom'),
        pl.col('POS').alias('pos'),
        ('SNP_' + pl.col('POS').cast(str) + '_' + pl.col('REF') + '_' +  pl.col('ALT')).alias('varname'),
        pl.lit(False).alias('is_STR'),
        pl.col('P').alias('p_val'),
        pl.col('BETA').alias('coeff'),
        pl.col('SE').alias('se'),
        pl.lit(df['method'][row]).alias('method'),
        pl.lit(df['region'][row]).alias('region'),
        pl.lit(df['replicate'][row]).alias('replicate').cast(int),
        pl.lit(vars_).alias('causal_vars'),
        pl.lit(betas).alias('causal_betas'),
        pl.lit(vars_betas).alias('causal_vars_and_betas'),
    ]))

    with open(df['str_assoc'][row]) as str_assoc_file:
        # confirm csv is not empty
        try:
            next(str_assoc_file) # header
            next(str_assoc_file) # potential first line of content
        except StopIteration:
            continue

    str_assoc_dfs.append(pl.scan_csv(
        df['str_assoc'][row],
        sep='\t'
    ).select([
        'chrom',
        'pos',
        ('STR_' + pl.col('pos').cast(str)).alias('varname'),
        pl.lit(True).alias('is_STR'),
        pl.col('p_simulation').alias('p_val'),
        pl.col('coeff_simulation').alias('coeff'),
        pl.col('se_simulation').alias('se'),
        pl.lit(df['method'][row]).alias('method'),
        pl.lit(df['region'][row]).alias('region'),
        pl.lit(df['replicate'][row]).alias('replicate').cast(int),
        pl.lit(vars_).alias('causal_vars'),
        pl.lit(betas).alias('causal_betas'),
        pl.lit(vars_betas).alias('causal_vars_and_betas'),
    ]))

assoc_df = pl.concat([*str_assoc_dfs, *snp_assoc_dfs])
finemap_df = pl.concat(finemap_dfs)
susie_df = pl.concat(susie_dfs)

final_df = susie_df.join(
    finemap_df,
    how='outer',
    on = ['varname', 'method', 'chrom', 'region', 'replicate', 'causal_vars', 'causal_betas', 'causal_vars_and_betas']
).join(
    assoc_df,
    how = 'left',
    on = ['varname', 'method', 'chrom', 'region', 'replicate', 'causal_vars', 'causal_betas', 'causal_vars_and_betas']
).collect()

final_df.write_csv('simulations_df.tab', sep='\t')
