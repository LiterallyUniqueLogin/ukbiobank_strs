#!/usr/bin/env python3

import argparse
import time

import glob
import os
import os.path

import bokeh.io
import bokeh.plotting
import bokeh.transform
import bokeh.util.hex
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

ukb = os.environ['UKB']

corr_cutoff = .8
p_val_thresh = 5e-8
pip_threshold = .8
# TODO finemap cutoff
# finemap_pip_cutoff = .8

def linear_int_interpolate(c1, c2, dist):
    c_new = []
    for coord1, coord2 in zip(c1, c2):
        c_new.append(coord1 + round((coord2 - coord1)*dist))
    return c_new

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotypes', nargs='+')
    phenotypes = parser.parse_args().phenotypes

    all_dfs = []
    susie_cs_min_abs_corrs = []
    finemap_cs_coverages = []
    unconverged_regions = []
    #underexplored_regions = []
    unfinished_regions = []

    for phenotype in phenotypes:
    
        pheno_dfs =  []
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

        assocs = pl.concat([str_assocs, snp_assocs]).filter(pl.col('p_val') <= p_val_thresh)

        regions_df = pl.read_csv(
            f'{ukb}/signals/regions/{phenotype}.tab',
            sep='\t'
        )
        for chrom, start, end, any_strs in zip(regions_df['chrom'], regions_df['start'], regions_df['end'], regions_df['any_strs']):
            if not any_strs:
                continue
            converged_fname = f'{ukb}/finemapping/susie_results/{phenotype}/{chrom}_{start}_{end}/converged.txt'
            if not os.path.exists(converged_fname):
                unfinished_regions.append((phenotype, chrom, start, end))
                continue
            with open(converged_fname) as converged_file:
                if not next(converged_file).strip() == 'TRUE':
                    unconverged_regions.append((phenotype, chrom, start, end))
                    continue
            print(f'Loading {phenotype} region {chrom}:{start}-{end}', flush=True)
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
                'susie_cs': [-1]*len(susie_vars),
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
            ).with_columns([
                pl.col('varname').str.extract('^[^_]*_([^_]*)', 1).cast(int).alias('pos'),
                pl.col('varname').str.extract('^[^_]*_[^_]*_([^_]*)_.*', 1).str.lengths().cast(int).alias('reflen'),
                pl.col('varname').str.extract('^[^_]*_[^_]*_[^_]*_([^_]*)', 1).str.lengths().cast(int).alias('altlen'),
                pl.col('varname').str.contains('^STR').alias('is_STR'),
                pl.lit(f'{phenotype}_{chrom}_{start}_{end}').alias('region'),
                pl.lit(chrom).alias('chrom').cast(int),
                pl.lit(phenotype).alias('phenotype')
            ]).sort('susie_idx')

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
                if min_abs_corr < corr_cutoff:
                    continue
                real_cs_count += 1
                # could worry about variants being in multiple CSes
                df = df.with_column(
                    pl.when(pl.col('susie_idx').is_in(cs_susie_idx))
                      .then(cs_id)
                      .otherwise(pl.col('susie_cs'))
                      .alias('susie_cs')
                )
            pheno_dfs.append(df)
            '''
            if real_cs_count >= 10:
                underexplored_regions.append((phenotype, chrom, start, end))
            '''
        pheno_dfs = [df.select(
            pl.col('*').exclude('^alpha.*$')
        ) for df in pheno_dfs]
        pheno_df = pl.concat(pheno_dfs).join(
            assocs,
            how='left',
            on=['phenotype', 'chrom', 'is_STR', 'pos', 'reflen', 'altlen']
        ).collect()
        all_dfs.append(pheno_df)

    del df, susie_df, finemap_df, assocs, pheno_dfs, pheno_df
    susie_cs_min_abs_corrs = np.array(susie_cs_min_abs_corrs)
    finemap_cs_coverages = np.array(finemap_cs_coverages)

    total_df = pl.concat(all_dfs)
    #total_assocs = pl.concat(all_assocs).filter(pl.col('p_val') <= p_val_thresh)

    ''''
    start_time = time.time()
    print('Gathering data ... ', flush=True)
    total_df = total_df.join(
        total_assocs,
        how='left',
        on=['phenotype', 'chrom', 'is_STR', 'pos', 'reflen', 'altlen']
    ).collect()
    print(f'Done. Time: {time.time() - start_time:.2}')
    '''

    total_df.filter(
        ~pl.col('p_val').is_null() & (pl.col('p_val') <= p_val_thresh)
    ).to_csv(f'{ukb}/post_finemapping/intermediate_results/gathered_data.tab', sep='\t')

    print('Any vars with null Ps?', total_df.select(pl.col('p_val').is_null().alias('null?')).select(pl.any('null?').alias('any_nulls'))['any_nulls'][0])
    print('n regions', total_df.select(pl.col('region').unique().count().alias('region_count'))['region_count'][0])

    cses_per_region = total_df.filter(
        pl.col('susie_cs') >= 0
    ).filter(
        ~pl.col('p_val').is_null()
    ).groupby([
        'susie_cs', 'region'
    ]).agg(
        pl.col('p_val').min().alias('min_p'),
    ).filter(
        pl.col('min_p') <= p_val_thresh
    ).groupby('region').agg(
        pl.col('region').count().alias('n_cses')
    ).to_dict(False)['n_cses']
    print(f'avg cses (total PIP >= .9, min_p_val of CS members <= {p_val_thresh}) per region {np.mean(cses_per_region)}, ({np.std(cses_per_region)})')

    for filter_, text in ((pl.lit(True), ''), (pl.col('is_STR'), ' STR'), (~pl.col('is_STR'), ' SNP')):
        susie_hits_per_region = total_df.filter(
            filter_
        ).with_column(
            ((pl.col('susie_cs') >= 0) & (pl.col('susie_pip') >= pip_threshold) & (pl.col('p_val') <= p_val_thresh)).alias('susie_hit')
        ).groupby('region').agg(
            pl.col('susie_hit').sum().alias('n_susie_hits')
        ).to_dict(False)['n_susie_hits']
        print(f'avg susie{text} hits (var is in a CS, PIP >= {pip_threshold}, p_val <= {p_val_thresh}) per region {np.mean(susie_hits_per_region)}, ({np.std(susie_hits_per_region)})')

        finemap_hits_per_region = total_df.filter(
            filter_
        ).with_column(
            ((pl.col('finemap_pip') >= pip_threshold) & (pl.col('p_val') <= p_val_thresh)).alias('finemap_hit')
        ).groupby('region').agg(
            pl.col('finemap_hit').sum().alias('n_finemap_hits')
        ).select('n_finemap_hits').to_numpy()
        print(f'avg finemap{text} hits (PIP >= {pip_threshold}, p_val <= {p_val_thresh}) per region {np.mean(finemap_hits_per_region)}, ({np.std(finemap_hits_per_region)})')

        print('Exporting FINEMAP vs SuSiE PIP plots', flush=True)
        comparison_thresh = 0.3
        title = f'{text} with p-val <= {p_val_thresh} where at least one of SuSiE or FINEMAP PIP >= {comparison_thresh}'
        if text == '':
            title = 'Vars ' + title
        fig = bokeh.plotting.figure(
            width=1200,
            height=1200,
            title=title,
            x_axis_label='FINEMAP PIPs',
            y_axis_label='SuSiE PIPs',
        )
        fig.title.text_font_size = '30px'
        fig.axis.axis_label_text_font_size = '26px'
        fig.axis.major_label_text_font_size = '20px'

        fig.background_fill_color = None
        fig.border_fill_color = None
        fig.ygrid.grid_line_color = None
        fig.xgrid.grid_line_color = None
        fig.toolbar.logo = None
        fig.toolbar_location = None
        print(total_df.filter(filter_))
        print(total_df.filter(filter_ & (pl.col('p_val') <= p_val_thresh)))
        pips = total_df.filter(
            filter_ &
            (pl.col('p_val') <= p_val_thresh) &
            ((pl.col('finemap_pip') >= comparison_thresh) | (
                (pl.col('susie_pip') >= comparison_thresh) & (pl.col('susie_cs') >= 0)
            ))
        ).select(['susie_pip', 'finemap_pip'])
        print(pips)

        bin_size = .05
        bins = bokeh.util.hex.hexbin(pips['finemap_pip'].to_numpy().reshape(-1), pips['susie_pip'].to_numpy().reshape(-1), size=bin_size)

        palette = [
            linear_int_interpolate((134,204,195), (9,41,46), i/254) for i in range(-1, 255)
        ]
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
            color_mapper = color_mapper,
            width=70,
            major_label_text_font_size = '20px'
        )
        fig.add_layout(color_bar, 'right')
        ext=text.replace(' ', '_')
        bokeh.io.export_png(fig, filename=f'{ukb}/export_scripts/results/finemap_pip_vs_susie_pip{ext}.png')
        bokeh.io.export_svg(fig, filename=f'{ukb}/export_scripts/results/finemap_pip_vs_susie_pip{ext}.svg')


    print(f'unconverged regions: {unconverged_regions}')
    print(f'unfinished regions: {unfinished_regions}')
    #print(f'underexplored regions: {underexplored_regions}')

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
    ys = [np.sum((left_edge <= susie_cs_min_abs_corrs) & (susie_cs_min_abs_corrs < left_edge + step)) for left_edge in left_edges]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+step)

    print('Exporting cs plots', flush=True)
    bokeh.io.export_png(fig, filename=f'{ukb}/export_scripts/results/cs_min_abs_corrs.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/export_scripts/results/cs_min_abs_corrs.svg')

    fig = bokeh.plotting.figure(
        width=1200,
        height=1200,
        title=f'Number of SuSie CSes min absolute corr >= {corr_cutoff} per region',
        x_axis_label='# cses in the region',
        y_axis_label='# regions',
    )
    fig.axis.axis_label_text_font_size = '30px'
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color = None
    fig.toolbar_location = None
    left_edges = np.arange(0, max(cses_per_region)+1)
    ys = [np.sum((left_edge <= cses_per_region) & (cses_per_region < left_edge + 1)) for left_edge in left_edges]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+1)

    print('Exporting cs per region plots', flush=True)
    bokeh.io.export_png(fig, filename=f'{ukb}/export_scripts/results/cses_per_region.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/export_scripts/results/cses_per_region.svg')

    fig = bokeh.plotting.figure(
        width=1200,
        height=1200,
        title=f'Number of FINEMAP vars with PIP >= {pip_threshold} per region',
        x_axis_label='# hits in the region',
        y_axis_label='# regions',
    )
    fig.axis.axis_label_text_font_size = '30px'
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color = None
    fig.toolbar_location = None
    left_edges = np.arange(0, max(finemap_hits_per_region)+1)
    ys = [np.sum((left_edge <= finemap_hits_per_region) & (finemap_hits_per_region < left_edge + 1)) for left_edge in left_edges]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+1)

    print('Exporting finemap hits per region plots', flush=True)
    bokeh.io.export_png(fig, filename=f'{ukb}/export_scripts/results/finemap_hits_per_region.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/export_scripts/results/finemap_hits_per_region.svg')

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
    include = susie_cs_min_abs_corrs >= corr_cutoff
    max_total_pip = max(1, np.max(finemap_cs_coverages[include]))
    step = 0.01
    left_edges = np.arange(0, max_total_pip + step, step)
    ys = [np.sum((left_edge <= finemap_cs_coverages[include]) & (finemap_cs_coverages[include] < left_edge + step)) for left_edge in left_edges]
    fig.quad(top=ys, bottom=0, left=left_edges, right=left_edges+step)

    print('Exporting FINEMAP CS PIP plots', flush=True)
    bokeh.io.export_png(fig, filename=f'{ukb}/export_scripts/results/susie_cs_finemap_total_pips.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/export_scripts/results/susie_cs_finemap_total_pips.svg')

    total_cses = np.sum(include)
    total_cses_large_finemap_pip = np.sum(finemap_cs_coverages[include] >= pip_threshold)
    print(f'SuSiE CSes with min_abs_corr >= {corr_cutoff} with FINEMAP total PIP >= {pip_threshold}: {total_cses_large_finemap_pip} ({total_cses_large_finemap_pip/total_cses:%})')

    susie_pip_threshold_for_finemap = .3
    n_replicates_from_finemap = total_df.filter(
        (pl.col('susie_cs') >= 0)  &
        (pl.col('susie_pip') >= susie_pip_threshold_for_finemap) &
        (pl.col('finemap_pip') >= pip_threshold)
    ).shape[0]
    n_finemap_total = total_df.filter(
        pl.col('finemap_pip') >= pip_threshold
    ).shape[0]
    print(f'FINEMAP hits with PIP >= {pip_threshold} in a SuSiE CS with abs corr >= {corr_cutoff} and SuSiE PIP >= {susie_pip_threshold_for_finemap}: {n_replicates_from_finemap} ({n_replicates_from_finemap/n_finemap_total:%})')

    for (curr_df,text) in [(total_df, 'all hits no filter'), (total_df.filter(pl.col('p_val') <= 1e-10), 'all hits p<=1e-10')]:
        print(text)
        var_thresh1 = .8
        var_thresh2 = .3
        for susie_thresh in (var_thresh1, var_thresh2):
            for finemap_thresh in (var_thresh1, var_thresh2):
                count = curr_df.filter(
                    (pl.col('susie_cs') >= 0) &
                    (pl.col('susie_pip') >= susie_thresh) &
                    (pl.col('finemap_pip') >= finemap_thresh)
                ).shape[0]
                print(f'Vars in a SuSiE CS with SuSIE PIP >= {susie_thresh} and with FINEMAP PIP >= {finemap_thresh}: {count}')
            
        for susie_thresh in (var_thresh1, var_thresh2):
            count = curr_df.filter(
                (pl.col('susie_cs') >= 0) &
                (pl.col('susie_pip') >= susie_thresh) &
                (pl.col('finemap_pip') < var_thresh2)
            ).shape[0]
            print(f'Vars in a SuSiE CS with SuSIE PIP >= {susie_thresh} with FINEMAP PIP < {var_thresh2}: {count}')
        for finemap_thresh in (var_thresh1, var_thresh2):
            count = curr_df.filter(
                (pl.col('finemap_pip') >= finemap_thresh) & (
                    (pl.col('susie_cs') < 0) |
                    (pl.col('susie_pip') < var_thresh2)
                )
            ).shape[0]
            print(f'Vars with FINEMAP PIP >= {finemap_thresh} either not in a SuSiE CS or having SuSiE PIP <= {var_thresh2}: {count}')

    # Not going to report susie alphas v pips - just know that they're similar if we look
    # at vars in good credible sets and not otherwise
    '''
    for vartype in '', 'STR', 'SNP':
        for all_ in True, False:
            # susie v finemap
            figure, ax = plt.subplots(1, 1, figsize=(12,12))

            if vartype == '':
                title_inset = 'all'
            else:
                title_inset = vartype
            if not all_:
                title_inset += f' credible-set-var (min_abs_corr >= {corr_cutoff})'
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
    main()
