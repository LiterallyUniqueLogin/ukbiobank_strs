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
import matplotlib.cm
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import upsetplot

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

def load_finemap(results_regions_dir, regions = None):
    unfinished_regions = []
    underexplored_regions = []
    dfs = []

    if regions is None:
        regions = []
        dirlist = os.listdir(results_regions_dir)
        for dir_ in dirlist:
            match = re.match('^([0-9]+)_[0-9]+_[0-9]+$', dir_)
            if not match:
                continue
            regions.append((dir_, match[1]))
    for (region, chrom) in regions:
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
    ]).distinct().filter(~(
        (
            (pl.col('phenotype') == 'total_bilirubin') &
            (pl.col('signal_region') == '12_19976272_22524428')
        ) |
        (
            (pl.col('phenotype') == 'mean_platelet_volume') &
            (pl.col('signal_region') == '17_2341352_2710113')
        )
    )).collect().to_dict(False)

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
    #for count, (phenotype, regions) in [(0, ('aspartate_aminotransferase', phenos_to_regions['aspartate_aminotransferase'])), (1, ('total_bilirubin', phenos_to_regions['total_bilirubin']))]:
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
        original_finemaps     = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}', regions=regions)
        print('ratio FINEMAP')
        ratio_finemaps        = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}.snp_str_ratio_4', regions=regions)
        print('total prob FINEMAP')
        total_prob_finemaps   = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}.total_prob_4', regions=regions)
        print('prior std FINEMAP')
        prior_std_finemaps    = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}.prior_std_0.005', regions=regions)
        print('conv tol FINEMAP')
        conv_tol_finemaps     = load_finemap(f'{ukb}/finemapping/finemap_results/{phenotype}.prob_conv_sss_tol_0.0001', regions=regions)
        print('mac FINEMAP')
        mac_finemaps          = load_finemap(f'{ukb}/finemapping/finemap_results_mac_100/{phenotype}', regions=regions)
        print('gt threshold FINEMAP')
        gt_threshold_finemaps = load_finemap(f'{ukb}/finemapping/finemap_results_threshold_0.0005/{phenotype}', regions=regions)

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
            prior_std_finemaps,
            how='outer',
            on=['chrom', 'varname'],
            suffix='_prior_std'
        ).drop('region_prior_std').join(
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
        'finemap_pip',
        'finemap_pip_ratio',
        'finemap_pip_total_prob',
        'finemap_pip_prior_std',
        'finemap_pip_conv_tol',
        'finemap_pip_mac',
        'finemap_pip_gt_thresh',
    ])
    total_df.to_csv(f'{ukb}/post_finemapping/intermediate_results/finemapping_concordance.tab', sep='\t')
    assert total_df.select((pl.col('region') == '').any().alias('region'))['region'].to_numpy()[0] == False
    assert np.all(~np.isnan(total_df['p_val'].to_numpy()))
    assert np.all(~np.isnan(total_df['finemap_pip'].to_numpy()))
    assert np.all(np.isnan(total_df['susie_alpha'].to_numpy()) == np.isnan(total_df['susie_alpha_hardcall'].to_numpy()))
    assert np.all(np.isnan(total_df['susie_alpha'].to_numpy()) == np.isnan(total_df['susie_alpha_ratio'].to_numpy()))
    assert np.all(1 == total_df.groupby(['phenotype', 'chrom', 'varname']).agg([pl.count()]).sort('count')['count'].to_numpy())
    return total_df

def linear_int_interpolate(c1, c2, dist):
    c_new = []
    for coord1, coord2 in zip(c1, c2):
        c_new.append(coord1 + round((coord2 - coord1)*dist))
    return c_new

def putatively_causal_hits_comparison(regenerate):
    if regenerate:
        total_df = putatively_causal_hits_df()
    else:
        total_df = pl.read_csv(f'{ukb}/post_finemapping/intermediate_results/finemapping_concordance.tab', sep='\t')

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
    ).drop('is_STR').to_csv(f'{ukb}/post_finemapping/intermediate_results/original_set.tab', sep='\t')

    susie_cols = total_df.select([
        pl.col('^susie_alpha.*$'),
    ]).columns
    finemap_cols = total_df.select([
        pl.col('^finemap_pip.*$')
    ]).columns
    print(susie_cols + finemap_cols)
    print(len(susie_cols + finemap_cols))
    print('n pass threshold .8 or more in original runs and each replicate with distinct loci:')
    pass_all_threshes = total_df.filter(
        pl.col('is_STR') &
        (pl.col('p_val') <= 1e-10) &
        (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols if 'ratio' not in col]) == 8)
        #(pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols]) == 10)
    )
    pass_all_threshes.select(['phenotype', 'region', 'chrom', 'pos', 'p_val']).to_csv(f'{ukb}/post_finemapping/intermediate_results/best_set.tab', sep='\t')
    print(pass_all_threshes.select(['pos', 'chrom']).shape[0])
    exit()

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
            plt.suptitle(f'Fine-mapping conditions breakdown')
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

    exit()
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
        bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/{mapper.lower()}_consistency_{suffix}.png')
        bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/{mapper.lower()}_consistency_{suffix}.svg')

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
