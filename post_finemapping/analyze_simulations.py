#!/usr/bin/env python3

import argparse

import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('simulations_df')
parser.add_argument('flank_start_to_start_and_end_pos')
args = parser.parse_args()

df = pl.read_csv(
    args.simulations_df,
    separator='\t',
    dtypes = {'causal_betas': str}
)
fix_hipstr  = pl.read_csv(
    args.flank_start_to_start_and_end_pos,
    separator='\t'
)

for colname in df.columns:
    if colname not in ['susie_alpha', 'susie_cs']:
        assert df.filter(pl.col(colname).is_null()).shape[0] == 0

df = df.with_columns([
    pl.col('varname') + '_',
    pl.col('causal_vars').str.replace_all('\t', '_\t') + '_',
    pl.when(~pl.col('susie_alpha').is_null()).then(pl.col('susie_alpha')).otherwise(pl.lit(0)),
    pl.when(~pl.col('susie_cs').is_null()).then(pl.col('susie_cs')).otherwise(pl.lit(-1)),
])

for cvars, chrom in df.unique(['method', 'chrom', 'region', 'replicate']).select(['causal_vars', 'chrom']).iter_rows():
    for cvar in cvars.split('\t'):
        if cvar[:3] == 'STR':
            fix = fix_hipstr.filter((pl.col('snpstr_pos') == int(cvar.split('_')[1])) & (pl.col('chrom') == chrom))
            assert fix.shape[0] in {0, 1}
            if fix.shape[0] == 1:
                #print(f'Swapping {cvar} {fix[0, "pos"]}', flush=True)
                df = df.with_columns(pl.when(pl.col('chrom') == chrom).then(pl.col('causal_vars').str.replace_all(cvar, f'STR_{fix[0, "pos"]}_')).otherwise(pl.col('causal_vars')))

assert df.unique(['method', 'chrom', 'region', 'replicate', 'varname']).shape[0] == df.shape[0]
assert df.unique('method').shape[0] == 5
assert df.unique('replicate').shape[0] == 3
assert df.select(['method', 'chrom', 'region', 'replicate']).unique().groupby(['method', 'chrom', 'region']).agg(pl.count()).select((pl.col('count') == 3).all()).item()
print('3 reps per region')

methods = ['random_one_var', 'random_two_var', 'random_three_var', 'susie_no_strs', 'susie_snps_strs']
assert set(methods) == set(df.unique('method')['method'])
#methods = ['susie_snps_strs']

stat_names = []
stats = []

def register_stat(stat_name, stat):
    stat_names.append(stat_name)
    stats.append(stat)

register_stat(
    'n regions',
    lambda sub: sub.unique(['region', 'chrom']).shape[0]
)
register_stat(
    'n simulations no c or non-c found',
    lambda sub: 3*sub.unique(['region', 'chrom']).shape[0] - sub.filter(pl.col('found')).select(['region', 'chrom', 'replicate']).unique().shape[0]
)
register_stat(
    'avg cvars per region',
    lambda sub: '{:.2f}'.format(
        sub.unique(['region', 'chrom']).select(pl.col('causal_vars').str.split('\t').list.lengths().mean()).item()
    )
)
register_stat(
    'avg cSTRs per region',
    lambda sub: '{:.2f}'.format(
        sub.unique(['region', 'chrom']).select(pl.col('causal_vars').str.split('\t').list.eval(pl.element().str.contains('STR').cast(int)).list.sum().mean()).item()
    ) if sub.select(pl.col('causal_vars').str.contains('STR').any()).item() else 'NA'
)
register_stat(
    'cSNPs found',
    lambda sub: sub.filter(pl.col('found') & ~pl.col('is_STR') & pl.col('causal_vars').str.contains(pl.col('varname'))).shape[0]
)
register_stat(
    'non-cSNPs found',
    lambda sub: sub.filter(pl.col('found') & ~pl.col('is_STR') & ~pl.col('causal_vars').str.contains(pl.col('varname'))).shape[0]
)
register_stat(
    '% cSNPs found',
    lambda sub: '{:.2f}'.format(
        sub.filter(pl.col('found') & ~pl.col('is_STR') & pl.col('causal_vars').str.contains(pl.col('varname'))).shape[0] / \
        sub.unique(['region', 'chrom', 'replicate']).select(pl.col('causal_vars').str.split('\t').list.eval(pl.element().str.contains('SNP').cast(int)).list.sum().sum()).item()
    ) if sub.select(pl.col('causal_vars').str.contains('SNP').any()).item() else 'NA'
)
register_stat(
    'cSTRs found',
    lambda sub: sub.filter(pl.col('found') & pl.col('is_STR') & pl.col('causal_vars').str.contains(pl.col('varname'))).shape[0] if sub.select(pl.col('causal_vars').str.contains('STR').any()).item() else 'NA'
)
register_stat(
    'non-cSTRs found',
    lambda sub: sub.filter(pl.col('found') & pl.col('is_STR') & ~pl.col('causal_vars').str.contains(pl.col('varname'))).shape[0]
)
register_stat(
    '% cSTRs found',
    lambda sub: '{:.2f}'.format(
        sub.filter(pl.col('found') & pl.col('is_STR') & pl.col('causal_vars').str.contains(pl.col('varname'))).shape[0] / \
        sub.unique(['region', 'chrom', 'replicate']).select(pl.col('causal_vars').str.split('\t').list.eval(pl.element().str.contains('STR').cast(int)).list.sum().sum()).item()
    ) if sub.select(pl.col('causal_vars').str.contains('STR').any()).item() else 'NA'
)

#print(
#    'method', 'cSNPs found', 'non-cSNPs misfound', '% cSNPs found', 'cSTRs found', 'non-cSTRs found', '% cSTRs found', 'n sims no cvars found', 'avg cvars per region', 'n regions', 'avg replicates per region', sep='\t'
#)

df = df.with_columns()

print('criteria', 'method', *stat_names, sep='\t')
for found_name, found_expr in [
    ('both CP >= 0.8', (pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 1) & (pl.col('finemap_pip') >= 0.8)),
    ('FINEMAP CP >= 0.8', pl.col('finemap_pip') >= 0.8),
    ('FINEMAP CP >= 0.5', pl.col('finemap_pip') >= 0.5),
    ('SuSiE CP >= 0.8', (pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 1)),
    ('SuSiE CP >= 0.5', (pl.col('susie_alpha') >= 0.5) & (pl.col('susie_cs') >= 1)),
    ('Best in SuSiE CS', (pl.col('susie_cs') >= 1) & (pl.col('susie_alpha') == pl.col('susie_alpha').max().over(['region', 'chrom', 'replicate', 'susie_cs']))),
#    ('In SuSiE CS', pl.col('susie_cs') >= 1)
]:
    for method in methods:
        found_df = df.filter(pl.col('method') == method).with_columns(found_expr.alias('found'))
        print(found_name, method, *[stat(found_df) for stat in stats], sep='\t')

for found_name, found_expr in [
    ('both CP >= 0.8', (pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 1) & (pl.col('finemap_pip') >= 0.8)),
    ('FINEMAP CP >= 0.8', pl.col('finemap_pip') >= 0.8),
    ('FINEMAP CP >= 0.5', pl.col('finemap_pip') >= 0.5),
    ('SuSiE CP >= 0.8', (pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 1)),
    ('SuSiE CP >= 0.5', (pl.col('susie_alpha') >= 0.5) & (pl.col('susie_cs') >= 1)),
    ('Best in SuSiE CS', (pl.col('susie_cs') >= 1) & (pl.col('susie_alpha') == pl.col('susie_alpha').max().over(['region', 'chrom', 'replicate', 'susie_cs']))),
#    ('In SuSiE CS', pl.col('susie_cs') >= 1)
]:
    for method in methods:
        found_df = df.filter(pl.col('method') == method).with_columns((found_expr & (pl.col('p_val') < 5e-8)).alias('found'))
        print(f'GWsig & {found_name}', method, *[stat(found_df) for stat in stats], sep='\t')

