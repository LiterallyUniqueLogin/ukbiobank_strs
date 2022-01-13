#!/usr/bin/env python3

import argparse

import bokeh.plotting
import bokeh.util.hex
import numpy as np
import polars as pl
import scipy.stats

def linear_int_interpolate(c1, c2, dist):
    c_new = []
    for coord1, coord2 in zip(c1, c2):
        c_new.append(coord1 + round((coord2 - coord1)*dist))
    return c_new

parser = argparse.ArgumentParser()
parser.add_argument('outdir')
parser.add_argument('results_table', help='cols: chrom, pos, FINEMAP_pcausal, SuSiE_CS_pcausal')
parser.add_argument('pos_to_snpstr_pos', help='cols: chrom, pos, snpstr_pos')
parser.add_argument('chrom_tables', nargs = '+',
                    help='In chromosome order. 4 cols: pos, chance of length confusionn, avg abs length confusion, normalized avg abs length confusion')
args = parser.parse_args()
outdir = args.outdir
results_fname = args.results_table
chrom_fnames = args.chrom_tables

loci = pl.concat([
    pl.scan_csv(
        chrom_fname,
        sep='\t'
    ).with_column(pl.lit(chrom_num+1).cast(int).alias('chrom'))
    for chrom_num, chrom_fname in enumerate(chrom_fnames)
]).collect()

pos_to_snpstr_pos = pl.scan_csv(
    args.pos_to_snpstr_pos,
    sep='\t'
).collect()

cols = ['normalized_avg_abs_length_confusion', 'chance_of_length_confusion']
#loci_cols = [col for col in loci.columns if col != 'pos' and col != 'chrom']

results = pl.scan_csv(
    results_fname,
    sep='\t',
    null_values='NA'
).with_column(
    pl.when(pl.col('SuSiE_CS_pcausal').is_null()).then(0).otherwise(pl.col('SuSiE_CS_pcausal')).alias('SuSiE_pcausal')
).with_column(
    (pl.col('FINEMAP_pcausal')-pl.col('SuSiE_pcausal')).alias('discrepancy')
).with_column(
    pl.col('discrepancy').abs().alias('abs_discrepancy')
).filter(
    (pl.col('FINEMAP_pcausal') >= .8) | (pl.col('SuSiE_pcausal') >= .8)
).select(['chrom', 'start_pos', 'discrepancy', 'abs_discrepancy']).collect()

df = results.join(
    pos_to_snpstr_pos,
    how='left',
    left_on=['chrom', 'start_pos'],
    right_on=['chrom', 'pos']
).join(
    loci,
    how='left',
    left_on = ['chrom', 'snpstr_pos'],
    right_on = ['chrom', 'pos']
)

assert df.filter(pl.col(cols[0]).is_null()).shape[0] == 0

unique_df = df.groupby(['chrom', 'start_pos']).agg([
    pl.col(col).first().keep_name() for col in cols
])

print('Plotting discrepancies ...', flush=True)
for col, x_axis_label in (('discrepancy', 'discrepancy (positive means FINEMAP likes this locus more)'), ('abs_discrepancy', 'abs_discrepancy')):
    figure = bokeh.plotting.figure(
        title=f'Distribution of {col}',
        x_axis_label=x_axis_label,
        y_axis_label='pdf',
        width=1200,
        height=1200,
        tools=''
    )
    figure.background_fill_color = None
    figure.border_fill_color = None
    figure.grid.grid_line_color = None
    figure.title.text_font_size='30px'
    figure.axis.axis_label_text_font_size='26px'
    figure.axis.major_label_text_font_size='20px'

    if col == 'discrepancy':
        xs=np.arange(-1, 1+.001, .001)
    else:
        xs=np.arange(0, 1+.001, .001)
    figure.varea(
        x=xs,
        y1=np.zeros(xs.shape),
        y2=scipy.stats.gaussian_kde(df[col].to_numpy())(xs),
    )

    bokeh.io.export_png(figure, filename=f'{outdir}/{col}.png')
    bokeh.io.export_svg(figure, filename=f'{outdir}/{col}.svg')

for col in cols:
    print(f'Plotting column {col} after subsetting to hits ...', flush=True)
    max_val = loci.select(pl.col(col).max()).to_numpy()
    min_val = loci.select(pl.col(col).min()).to_numpy()
    n_steps = 1000
    step_size = (max_val - min_val)/n_steps
    xs = np.arange(min_val, max_val + step_size, step_size)
    ys_hits = scipy.stats.gaussian_kde(unique_df[col].to_numpy())(xs)
    ys_all = scipy.stats.gaussian_kde(loci[col].to_numpy())(xs)

    if col.startswith('chance'):
        unit = '%'
    elif col.startswith('avg'):
        unit = 'repeat units'
    else:
        assert col.startswith('normalized')
        unit = 'length standard deviations'
    x_axis_label = f'{col} ({unit})'
   
    ''' 
    figure = bokeh.plotting.figure(
        title=col.capitalize() + ' distribution after subsetting to hits',
        x_axis_label=x_axis_label,
        y_axis_label='pdf',
        width=1200,
        height=1200,
        tools=''
    )
    figure.background_fill_color = None
    figure.border_fill_color = None
    figure.grid.grid_line_color = None
    figure.title.text_font_size='30px'
    figure.axis.axis_label_text_font_size='26px'
    figure.axis.major_label_text_font_size='20px'

    figure.varea(
        x=xs,
        y1=np.zeros(xs.shape),
        y2=ys_hits,
        legend_label='hits',
        alpha=0.5,
        color='blue'
    )
    figure.varea(
        x=xs,
        y1=np.zeros(xs.shape),
        y2=ys_all,
        legend_label='all_imputed_loci',
        alpha=0.5,
        color='red'
    )

    figure.legend.label_text_font_size = '22px'

    bokeh.io.export_png(figure, filename=f'{outdir}/{col}_distribution_for_hits.png')
    bokeh.io.export_svg(figure, filename=f'{outdir}/{col}_distribution_for_hits.svg')
    '''

    figure = bokeh.plotting.figure(
            title=col.capitalize() + ' vs abs finemapping discrepency',
            x_axis_label=x_axis_label,
            y_axis_label='abs finemapping discrepency',
            width=1200,
            height=1200,
            tools=''
    )
    figure.background_fill_color = None
    figure.border_fill_color = None
    figure.grid.grid_line_color = None
    figure.title.text_font_size='30px'
    figure.axis.axis_label_text_font_size='26px'
    figure.axis.major_label_text_font_size='20px'


    max_x_val = df.select(pl.col(col).max()).to_numpy()
    min_x_val = df.select(pl.col(col).min()).to_numpy()
    max_y_val = df.select(pl.col('abs_discrepancy').max()).to_numpy()
    min_y_val = df.select(pl.col('abs_discrepancy').min()).to_numpy()
    aspect_ratio = float((max_y_val-min_y_val)/(max_x_val - min_x_val))

    bin_size=.05
    bins = bokeh.util.hex.hexbin(
        df[col].to_numpy().reshape(-1),
        df['abs_discrepancy'].to_numpy().reshape(-1),
        size=bin_size,
        aspect_scale = aspect_ratio
    )
    palette = [
        linear_int_interpolate((134,204,195), (9,41,46), i/254) for i in range(-1, 255)
    ]
    cmap = bokeh.transform.linear_cmap(
        'counts',
        palette = palette,
        low=1,
        high=max(bins.counts),
        low_color=(255, 255, 255)
    )
    color_mapper = bokeh.models.LinearColorMapper(
        palette = palette,
        low=1,
        high=max(bins.counts)
    )

    figure.hex_tile(q='q', r='r', size=bin_size, line_color=None, source=bins, fill_color=cmap)
    color_bar = bokeh.models.ColorBar(
        color_mapper = color_mapper,
        width=70,
        major_label_text_font_size = '20px'
    )
    figure.add_layout(color_bar, 'right')

    bokeh.io.export_png(figure, filename=f'{outdir}/{col}_vs_finemapping_discrepancy.png')
    bokeh.io.export_svg(figure, filename=f'{outdir}/{col}_vs_finemapping_discrepancy.svg')

print('Done')
