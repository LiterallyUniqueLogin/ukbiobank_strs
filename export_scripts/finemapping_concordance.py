#!/usr/bin/env python3

import argparse
import glob
import os
import os.path

import bokeh.io
import bokeh.plotting
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

ukb = os.environ['UKB']

cutoff = .8

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotypes', nargs='+')
    phenotypes = parser.parse_args().phenotypes

    all_dfs = []
    susie_cs_min_abs_corrs = []
    finemap_cs_coverages = []
    unconverged_regions = []
    underexplored_regions = []

    n_regions = 0
    for phenotype in phenotypes:
        regions_df = pl.read_csv(
            f'{ukb}/signals/regions/{phenotype}.tab',
            sep='\t'
        )
        for chrom, start, end, any_strs in zip(regions_df['chrom'], regions_df['start'], regions_df['end'], regions_df['any_strs']):
            if not any_strs:
                continue
            with open(f'{ukb}/finemapping/susie_results/{phenotype}/{chrom}_{start}_{end}/converged.txt') as converged_file:
                if not next(converged_file).strip() == 'TRUE':
                    unconverged_regions.append((phenotype, chrom, start, end))
                    continue
            print(f'Loading {phenotype} region {chrom}:{start}-{end}', flush=True)
            n_regions += 1
            with open(f'{ukb}/finemapping/susie_results/{phenotype}/{chrom}_{start}_{end}/colnames.txt') as var_file:
                susie_vars = [line.strip() for line in var_file]
            alphas = pl.scan_csv(
                f'{ukb}/finemapping/susie_results/{phenotype}/{chrom}_{start}_{end}/alpha.tab',
                sep='\t',
                has_header=False
            ).collect().to_numpy().T
            n_alphas = alphas.shape[1]
            susie_pips=1-np.prod(1-alphas, axis=1)
            assert susie_pips.shape[0] == len(susie_vars)
            susie_idx = np.arange(len(susie_vars)) + 1
            susie_df = pl.DataFrame({
                'varname': susie_vars,
                'susie_pip': susie_pips,
                'susie_alpha': np.zeros(len(susie_vars)),
                'susie_cs': [False]*len(susie_vars),
                'susie_idx': susie_idx,
                **{ f'alpha_{i}': alphas[:, i] for i in range(n_alphas) }
            }).lazy()
            finemap_df = pl.scan_csv(
                f'{ukb}/finemapping/finemap_results/{phenotype}/{chrom}_{start}_{end}/finemap_output.snp',
                sep=' '
            ).select([
                pl.col('rsid').alias('varname'),
                pl.col('prob').alias('finemap_pip')
            ])

            df = susie_df.join(
                finemap_df,
                how='inner',
                on=['varname']
            ).sort('susie_idx')

            real_cs_count = 0
            for cs_fname in glob.glob(f'{ukb}/finemapping/susie_results/{phenotype}/{chrom}_{start}_{end}/cs*.txt'):
                cs_id = int(cs_fname.split('cs')[-1].split('.')[0])
                with open(cs_fname) as cs_file:
                    # susie uses 1 based indexing, python uses 0
                    # make sure cs idxs are in increasing order
                    cs_susie_idx = np.array([int(idx) for idx in next(cs_file).strip().split()])
                    assert np.all(cs_susie_idx[1:] - cs_susie_idx[:-1] > 0)
                    cs_susie_idx = pl.Series('cs_susie_idx', cs_susie_idx)
                    next(cs_file) # skip cs credibility
                    min_abs_corr, _, _ = [float(idx) for idx in next(cs_file).strip().split()]
                susie_cs_min_abs_corrs.append(min_abs_corr)
                finemap_cs_coverages.append(
                    df.filter(pl.col('susie_idx').is_in(cs_susie_idx)).select(pl.col('finemap_pip').sum()).collect()
                )
                df = df.with_column(
                    pl.when(pl.col('susie_idx').is_in(cs_susie_idx))
                      .then(pl.when(pl.col(f'alpha_{cs_id-1}') > pl.col('susie_alpha'))
                              .then(pl.col(f'alpha_{cs_id-1}'))
                              .otherwise(pl.col('susie_alpha')))
                      .otherwise(pl.col('susie_alpha'))
                      .alias('susie_alpha')
                )
                if min_abs_corr < cutoff:
                    continue
                real_cs_count += 1
                df = df.with_column(
                    pl.when(pl.col('susie_idx').is_in(cs_susie_idx))
                      .then(True)
                      .otherwise(pl.col('susie_cs'))
                      .alias('susie_cs')
                )
            if real_cs_count >= 10:
                underexplored_regions.append((phenotype, chrom, start, end))
            all_dfs.append(df)
    del df, susie_df, finemap_df
    susie_cs_min_abs_corrs = np.array(susie_cs_min_abs_corrs)
    finemap_cs_coverages = np.array(finemap_cs_coverages)

    all_dfs = [df.select(
        pl.col('*').exclude('^alpha.*$')
    ) for df in all_dfs]
    total_df = pl.concat(all_dfs).collect()

    print(f'n regions {n_regions}')
    print(f'n CSes over cutoff {cutoff}: {sum(susie_cs_min_abs_corrs >= cutoff)}')
    print(f'n vars: {total_df.select(pl.col("varname").count())[-1].to_numpy()}')
    print(f'n cutoff cs vars: {total_df.select(pl.col("susie_cs").sum()).to_numpy()}')
    print(f'unconverged regions: {unconverged_regions}')
    print(f'underexplored regions: {underexplored_regions}')

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
    ys = [np.sum((left_edge < susie_cs_min_abs_corrs) & (susie_cs_min_abs_corrs < left_edge + step)) for left_edge in left_edges]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+0.01)

    print('Exporting cs plots', flush=True)
    bokeh.io.export_png(fig, filename=f'{ukb}/export_scripts/results/cs_min_abs_corrs.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/export_scripts/results/cs_min_abs_corrs.svg')

    fig = bokeh.plotting.figure(
        width=1200,
        height=1200,
        title=f'FINEMAP total PIPs for SuSiE CSes with min_abs_corr >= {cutoff}',
        x_axis_label='FINEMAP PIPs',
        y_axis_label='# credible sets',
    )
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.ygrid.grid_line_color = None
    fig.xgrid.grid_line_color = None
    fig.toolbar.logo = None
    fig.toolbar_location = None
    include = susie_cs_min_abs_corrs >= cutoff
    max_total_pip = max(1, np.max(finemap_cs_coverages[include]))
    step = 0.01
    left_edges = np.arange(0, max_total_pip + step, step)
    ys = [np.sum((left_edge < finemap_cs_coverages[include]) & (finemap_cs_coverages[include] < left_edge + step)) for left_edge in left_edges]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+0.01)

    print('Exporting finemap CS PIP plots', flush=True)
    bokeh.io.export_png(fig, filename=f'{ukb}/export_scripts/results/susie_cs_finemap_total_pips.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/export_scripts/results/susie_cs_finemap_total_pips.svg')

    for vartype in '', 'STR', 'SNP':
        for all_ in True, False:
            # susie v finemap
            figure, ax = plt.subplots(1, 1, figsize=(12,12))

            if vartype == '':
                title_inset = 'all'
            else:
                title_inset = vartype
            if not all_:
                title_inset += f' credible-set-var (min_abs_corr >= {cutoff})'
            ax.set_title(f'SuSiE vs FINEMAP {title_inset} PIPs')
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlabel('FINEMAP PIPs')
            ax.set_ylabel('SuSiE PIPs')
            idx = total_df['varname'].str.contains(f'^{vartype}')
            idx = idx & (total_df['susie_cs'] | all_)
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list(name='meh', colors=['#e5f5f9', '#99d8c9', '#2ca25f'])
            hexbin = ax.hexbin(total_df[idx, 'finemap_pip'].to_numpy().flatten(), total_df[idx, 'susie_pip'].to_numpy().flatten(), gridsize=100, cmap=cmap, bins='log')
            cb = figure.colorbar(hexbin, ax=ax)
            cb.set_label('counts')

            if vartype != '':
                fname_suffix = '_' + vartype
            else:
                fname_suffix = ''
            if not all_:
                fname_suffix += '_credible_set_vars'
            out_fname = f'{ukb}/export_scripts/results/finemap_v_susie{fname_suffix}'
            print(f'Exporting var plot {out_fname}.svg', flush=True)
            plt.savefig(f'{out_fname}.svg')
            print(f'Exporting var plot {out_fname}.png', flush=True)
            plt.savefig(f'{out_fname}.png')
       
            # susie pip vs alpha
            figure, ax = plt.subplots(1, 1, figsize=(12,12))

            if vartype == '':
                title_inset = 'all'
            else:
                title_inset = vartype
            if not all_:
                title_inset += f' credible-set-var (min_abs_corr >= {cutoff})'
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



if __name__ == '__main__':
    main()
