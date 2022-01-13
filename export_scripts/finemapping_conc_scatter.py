#!/usr/bin/env python3

import numpy as np
import polars as pl
import bokeh.models
import bokeh.plotting
import bokeh.transform
import scipy.stats

size = 600
def linear_int_interpolate(c1, c2, dist):
    c_new = []
    for coord1, coord2 in zip(c1, c2):
        c_new.append(coord1 + round((coord2 - coord1)*dist))
    return c_new

palette = [
    linear_int_interpolate((134,204,195), (9,41,46), i/254) for i in range(-1, 255)
]

def heat_map(fig, x_vals, y_vals, fname, y_min=0):
    step = 0.02
    x_vals = x_vals.to_numpy()
    y_vals = y_vals.to_numpy()
    xs = np.arange(0, 1 + step, step)
    ys = np.arange(y_min, 1 + step, step)
    xs, ys = np.meshgrid(xs, ys)
    xs = xs.reshape(-1)
    ys = ys.reshape(-1)
    counts = [np.sum((x <= x_vals) & (x_vals < x + step) & (y <= y_vals) & (y_vals <= y + step)) for (x,y) in zip(xs, ys)]
    max_count = np.max(counts)

    color_mapper = bokeh.models.LinearColorMapper(palette=palette,low=1,high=max_count,low_color=(255,255,255))
    color_bar = bokeh.models.ColorBar(color_mapper=color_mapper, width=30)
    cmap = bokeh.transform.linear_cmap('counts', palette=palette,low=1,high=max_count,low_color=(255,255,255))
    cds = bokeh.models.ColumnDataSource(dict(
        left=xs,
        right=xs+step,
        top=ys+step,
        bottom=ys,
        counts=counts
    ))
    fig.quad(left='left', right='right', top='top', bottom='bottom', color=cmap, source=cds)
    cb_title = bokeh.models.Title(text='number of loci', align='center')
    fig.add_layout(color_bar, 'right')
    fig.add_layout(cb_title, 'right')
    bokeh.io.export_png(fig, filename=fname)

def weighted_heat_map(fig, x_vals, y_vals, weights, colorbar_text, fname):
    step=0.02
    x_vals = x_vals.to_numpy()
    y_vals = y_vals.to_numpy()
    weights=weights.to_numpy()
    xs = np.arange(0, 1+step, step)
    ys = np.arange(0, 1+step, step)
    xs, ys = np.meshgrid(xs, ys)
    xs = xs.reshape(-1)
    ys = ys.reshape(-1)
    weights = np.array([np.mean(weights[np.where((x <= x_vals) & (x_vals < x + step) & (y <= y_vals) & (y_vals <= y + step))]) for (x,y) in zip(xs, ys)])
    weights[np.isnan(weights)] = -1

    color_mapper = bokeh.models.LinearColorMapper(palette=palette,low=0,high=np.max(weights))
    color_bar = bokeh.models.ColorBar(color_mapper=color_mapper, width=30)
    cmap = bokeh.transform.linear_cmap('weights', palette=palette,low=0,high=np.max(weights), low_color=(255,255,255))
    cds = bokeh.models.ColumnDataSource(dict(
        left=xs,
        right=xs+step,
        top=ys+step,
        bottom=ys,
        weights=weights
    ))
    fig.quad(left='left', right='right', top='top', bottom='bottom', color=cmap, source=cds)
    cb_title = bokeh.models.Title(text=colorbar_text, align='center')
    fig.add_layout(color_bar, 'right')
    fig.add_layout(cb_title, 'right')
    bokeh.io.export_png(fig, filename=fname)

def count_alleles(thresh):
    def func(series_list):
        allele_dist_str = series_list[0]
        counts = []
        for str_ in allele_dist_str:
            split = str_.split(' ')
            count = 0
            for freq in split[1::2]:
                if float(freq[:-1]) >= thresh:
                    count += 1
            counts.append(count)
        return pl.Series(counts)
    return func

def main():
    df = pl.scan_csv('post_finemapping/intermediate_results/gathered_data.tab', sep='\t').filter((pl.col('susie_pip') >= 0.3) | (pl.col('finemap_pip') >= 0.3))
    df = df.with_column(
        (pl.col('susie_pip') - pl.col('finemap_pip')).alias('susie_f_pip_diff')
    ).with_column(
        pl.col('susie_f_pip_diff').abs().alias('abs_pip_diff')
    )
    locus_summary_df = pl.concat([
        pl.scan_csv(
            f'export_scripts/intermediate_results/chr{chrom}_loci_summary.tab',
            sep='\t'
        )
        for chrom in range(1, 23)
    ]).select(['chr', 'pos', 'multiallelicness', 'allele_dist'])
    allele_threshes = (0.0004, 0.002, 0.01, 0.05)
    #allele_threshes = [0.01]
    df = df.join(
        locus_summary_df,
        how='left',
        #left_on=['chrom', 'snpstr_pos'],
        left_on=['chrom', 'pos'],
        right_on=['chr', 'pos']
    ).collect()

    snp_df = df.filter(~pl.col('is_STR'))
    str_df = df.filter(pl.col('is_STR'))
    assert not str_df.select(pl.col('multiallelicness').is_null().any()).to_numpy()[0]

    str_df = str_df.with_columns([
        pl.apply('allele_dist', count_alleles(thresh), pl.UInt32).alias(f'alleles_{thresh}')
        for thresh in allele_threshes
    ])
    confusions = pl.concat([pl.scan_csv(f'side_analyses/length_confusion/chr{i}.tab', sep='\t').with_column(pl.lit(i).alias('chrom').cast(int)) for i in range(1,23)]).collect()
    merged_df = str_df.join(confusions, how='left',on=['chrom', 'pos'])

    step = 0.05
    fig = bokeh.plotting.figure(title='STR PIP histogram', width=size, height=size, x_axis_label='PIP', y_axis_label = 'density', tools='', toolbar_location=None)
    xs = np.arange(0,1+step,step)
    fig.line(
        x=xs[:-1],
        #y=scipy.stats.gaussian_kde(arr)(xs),
        y=np.histogram(str_df['susie_pip'], bins=xs, density=True)[0],
        color='red',
        legend_label='SuSiE STRs'
    )
    fig.line(
        x=xs[:-1],
        #y=scipy.stats.gaussian_kde(arr)(xs),
        y=np.histogram(str_df['finemap_pip'], bins=xs, density=True)[0],
        color='blue',
        legend_label='FINEMAP STRs'
    )
    fig.line(
        x=xs[:-1],
        #y=scipy.stats.gaussian_kde(arr)(xs),
        y=np.histogram(snp_df['susie_pip'], bins=xs, density=True)[0],
        color='green',
        legend_label='SuSiE SNPs'
    )
    fig.line(
        x=xs[:-1],
        #y=scipy.stats.gaussian_kde(arr)(xs),
        y=np.histogram(snp_df['finemap_pip'], bins=xs, density=True)[0],
        color='purple',
        legend_label='FINEMAP SNPs'
    )
    bokeh.io.export_png(fig, filename='post_finemapping/results/pip_histogram.png')

    fig = bokeh.plotting.figure(title='STR PIP scatterplot', width=size, height=size, x_axis_label='FINEMAP PIP', y_axis_label = 'SuSiE PIP', tools='', toolbar_location=None)
    fig.circle(str_df['susie_pip'], str_df['finemap_pip'])
    bokeh.io.export_png(fig, filename='post_finemapping/results/str_comp_pip_scatter.png')

    fig = bokeh.plotting.figure(title='STR PIP heatmap', width=size, height=size, x_axis_label='FINEMAP PIP', y_axis_label = 'SuSiE PIP', match_aspect=True, tools='', toolbar_location=None)
    heat_map(fig, str_df['finemap_pip'], str_df['susie_pip'], 'post_finemapping/results/str_comp_pip_heatmap.png')

    fig = bokeh.plotting.figure(title='STR PIPs', width=size, height=size, x_axis_label='FINEMAP PIP', y_axis_label = 'SuSiE PIP', match_aspect=True, tools='', toolbar_location=None)
    weighted_heat_map(fig, merged_df['finemap_pip'], merged_df['susie_pip'], merged_df['chance_of_length_confusion'], 'average chance of misgenotyping per sample at any such locus', 'post_finemapping/results/str_comp_pip_chance_map.png')
    
    fig = bokeh.plotting.figure(title='STR PIPs', width=size, height=size, x_axis_label='FINEMAP PIP', y_axis_label = 'SuSiE PIP', match_aspect=True, tools='', toolbar_location=None)
    weighted_heat_map(fig, merged_df['finemap_pip'], merged_df['susie_pip'], merged_df['normalized_avg_abs_length_confusion'], 'average number of standard deviations of misgenotyping per sample at any such locus', 'post_finemapping/results/str_comp_pip_sd_map.png')

    fig = bokeh.plotting.figure(title='SNP PIP scatterplot', width=size, height=size, x_axis_label='FINEMAP PIP', y_axis_label = 'SuSiE PIP', tools='', toolbar_location=None)
    fig.circle(snp_df['susie_pip'], snp_df['finemap_pip'])
    bokeh.io.export_png(fig, filename='post_finemapping/results/snp_comp_pip_scatter.png')

    fig = bokeh.plotting.figure(title='SNP PIP heatmap', width=size, height=size, x_axis_label='FINEMAP PIP', y_axis_label = 'SuSiE PIP', match_aspect=True, tools='', toolbar_location=None)
    heat_map(fig, snp_df['finemap_pip'], snp_df['susie_pip'], 'post_finemapping/results/snp_comp_pip_heatmap.png')

    color_mapper = bokeh.models.LinearColorMapper(palette=palette,low=0,high=1)
    color_bar = bokeh.models.ColorBar(color_mapper=color_mapper, width=30)
    cmap = bokeh.transform.linear_cmap('foo', palette=palette,low=0,high=1)

    fig = bokeh.plotting.figure(title='STR PIP scatterplot', width=size, height=size, x_axis_label='FINEMAP PIP', y_axis_label = 'SuSiE PIP', tools='', match_aspect=True, toolbar_location=None)
    cb_title = bokeh.models.Title(text='chance a genotype call at this locus is wrong', align='center')
    fig.add_layout(color_bar, 'right')
    fig.add_layout(cb_title, 'right')
    cds = bokeh.models.ColumnDataSource(dict(x=merged_df['finemap_pip'], y=merged_df['susie_pip'], color=[linear_int_interpolate((134,204,195), (9,41,46), val) for val in merged_df['chance_of_length_confusion']]))
    fig.circle(x='x', y='y', color='color', source=cds)
    bokeh.io.export_png(fig, filename='post_finemapping/results/colored_str_comp_pip_scatter.png')

    step=0.05
    for thresh in allele_threshes:
        for pip_thresh in (0.3, 0.8):
            for xs, x_label, out_loc, title, col in [
                (np.arange(-1,1+step,step),
                 'SuSiE PIP - FINEMAP PIP',
                 f'post_finemapping/results/pip_diff_density_allele_thresh_{thresh}_pip_thresh_{pip_thresh}.png',
                 f'PIP diff, STR allele penetrance threshold = {thresh:.4}',
                 'susie_f_pip_diff',
                ),
                (np.arange(0,1+step,step),
                 'absolute PIP difference',
                 f'post_finemapping/results/pip_abs_diff_density_allele_thresh_{thresh}_pip_thresh_{pip_thresh}.png',
                 f'absolute PIP diff, STR allele penetrance threshold = {thresh:.4}',
                 'abs_pip_diff'
                )
            ]:
                filter_exp = (pl.col('susie_pip') >= pip_thresh)  | (pl.col('finemap_pip') >= pip_thresh)
                fig = bokeh.plotting.figure(title=title, width=size, height=size, x_axis_label=x_label, y_axis_label = 'density', tools='', toolbar_location=None)
                fig.line(
                    x=xs[:-1],
                    y=np.histogram(snp_df.filter(filter_exp)[col].to_numpy(), bins=xs, density=True)[0],
                    #y=scipy.stats.gaussian_kde(snp_df['susie_f_pip_diff'].to_numpy())(xs),
                    color='black',
                    legend_label=f'SNPs (n={snp_df.shape[0]})'
                )
                for count, color in ((2, 'brown'), (3, 'red'), (4, 'orange')):
                    arr = str_df.filter(filter_exp).filter(pl.col(f'alleles_{thresh}') == count)[col].to_numpy()
                    fig.line(
                        x=xs[:-1],
                        #y=scipy.stats.gaussian_kde(arr)(xs),
                        y=np.histogram(arr, bins=xs, density=True)[0],
                        color=color,
                        legend_label=f'{count}-allele STRs (n={arr.shape[0]})'
                    )
                arr = str_df.filter(filter_exp).filter(pl.col(f'alleles_{thresh}') >= 5)[col].to_numpy()
                fig.line(
                    x=xs[:-1],
                    #y=scipy.stats.gaussian_kde(arr)(xs),
                    y=np.histogram(arr, bins=xs, density=True)[0],
                    color='gold',
                    legend_label=f'STRs with at least 5 alleles (n={arr.shape[0]})'
                )
                fig.add_layout(bokeh.models.Title(
                    text=f'Variants with PIP at least {pip_thresh} for SuSiE or FINEMAP'
                ), 'below')
                bokeh.io.export_png(fig, filename=out_loc)

    fig = bokeh.plotting.figure(title='STR PIP diff', width=size, height=size, x_axis_label='multiallelicness', y_axis_label = 'SuSiE PIP - FINEMAP PIP', tools='', toolbar_location=None)
    heat_map(fig, str_df['multiallelicness'], str_df['susie_f_pip_diff'], 'post_finemapping/results/str_pip_diff_heatmap.png', y_min=-1)
    fig = bokeh.plotting.figure(title='STR PIP abs diff', width=size, height=size, x_axis_label='multiallelicness', y_axis_label = 'absolute PIP difference', tools='', toolbar_location=None)
    heat_map(fig, str_df['multiallelicness'], str_df['abs_pip_diff'], 'post_finemapping/results/str_pip_abs_diff_heatmap.png')

    fig = bokeh.plotting.figure(title='PIP abs diff', width=size, height=size, x_axis_label='multiallelicness', y_axis_label = 'absolute PIP difference', tools='', toolbar_location=None)

if __name__ == '__main__':
    main()
