#!/usr/bin/env python3

import argparse

import bokeh.io
import bokeh.models
import bokeh.models.tickers
import bokeh.plotting
import numpy as np
import polars as pl

import region_plots
import phenotypes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile')
    parser.add_argument('ext')
    #cols chrom, start_pos, phenotype, association_p_value, multiallelicness, repeat_unit, pcausal,
    #relation_to_gene, other_ethnic_effect_directions
    parser.add_argument('hits_table')
    parser.add_argument('snpstr_correspondence_table') #cols chrom, pos, end_pos, snpstr_pos
    #parser.add_argument('per_pheno_summary_tables', help='json dict from pheno to summary table')
    #parser.add_argument('per_pheno_results', help='json dict from pheno to results file')
    #parser.add_argument('units', help='json dict from pheno to unit')

    args = parser.parse_args()
    #summary_table_fnames = json.loads(args.per_pheno_summary_tables)
    # units = json.loads(args.units)

    # make sure same phenos everywhere
    #assert set(units) == set(per_pheno_results)

    df = pl.scan_csv(
        args.hits_table,
        sep='\t',
        dtype={'other_ethnic_association_ps': str}
    )
    snpstr_correspondence = pl.scan_csv(
        args.snpstr_correspondence_table,
        sep='\t'
    )
    df = df.join(snpstr_correspondence, how='left', left_on=['chrom', 'start_pos'], right_on=['chrom', 'pos'])
    df = df.rename({'start_pos': 'pos', 'chrom': 'chr'})

    df = df.with_columns([
        pl.when(pl.col('association_p_value') == 0)
            .then(300)
            .otherwise(-pl.col('association_p_value').log10())
            .alias('association_p_value'),
        pl.when(pl.col('phenotype').is_in(list(phenotypes.haematological_phenotypes)))
            .then('red')
            .otherwise('blue')
            .alias('color'),
        pl.when(pl.col('multiallelicness')*30 > 5)
            .then(pl.col('multiallelicness')*30)
            .otherwise(5)
            .alias('size'),
        pl.col('chr').count().over(['chr', 'pos']).alias('num_assocs'),
        #(pl.col('association_p_value').max().over(['chr', 'pos']) + 4).alias('box_max'),
        #(pl.col('association_p_value').min().over(['chr', 'pos']) - 4).alias('box_min'),
    ]).with_columns([
        (pl.col('size') + 8).alias('box_width'),
        (pl.col('association_p_value').max().over(['chr', 'pos']) + 5 + pl.col('size')).alias('box_top'),
        pl.when(pl.col('association_p_value').min().over(['chr', 'pos']) - 5 - pl.col('size') < 0)
         .then(0)
         .otherwise(pl.col('association_p_value').min().over(['chr', 'pos']) - 5 - pl.col('size')).alias('box_bottom')
    ]).with_columns([
        ((pl.col('box_bottom') + pl.col('box_top'))/2).alias('box_center'),
        (pl.col('box_top') - pl.col('box_bottom')).alias('box_height'),
    ])
    df = df.collect()
    # chr, pos, end_pos, snpstr_pos (start)
 
    # any column will do
    '''
    df['unit'] = 'NA'
    df['unit'] = df['unit'].cast(pl.Categorical)

    for phenotype, unit in units.items():
        df[df['phenotype'] == phenotype, 'unit'] = args.units[phenotype]
    '''

    results_plot = bokeh.plotting.figure(
        width=1200,
        height=400,
        title='Putatively causal STRs by position',
        x_axis_label='Chromosome',
        y_axis_label='-log10(p-value)',
        y_range = (0, 320),
        tools='reset, save'
    )

    results_plot.xgrid.ticker = bokeh.models.tickers.FixedTicker(
        ticks = np.cumsum(region_plots.chr_lens[:-1])
    )
    results_plot.ygrid.grid_line_color = None

    if args.ext == 'html':
        region_plots.add_x_navigate(results_plot)

    region_plots.full_genome_x_axis(results_plot)
    df = region_plots.full_genome_polars_df(df)

    line_source = bokeh.models.ColumnDataSource(dict(
        x=[0, sum(region_plots.chr_lens)],
        y=[0, 0]
    ))
    multi_assocs = df.filter(pl.col('num_assocs') > 1)

    results_plot.rect(
        x='plot_pos',
        y='box_center',
        height='box_height',
        width = 'box_width',
        width_units = 'screen',
        color='grey',
        #fill_color=None,
        source=bokeh.models.ColumnDataSource({col: multi_assocs[col] for col in multi_assocs.columns}),
        legend_label = 'Same signal multiple phenos'
    )

    # x axis
    results_plot.line(
        x='x',
        y='y',
        source=line_source,
        color='black',
        line_width=1
    )


    '''
    results_plot.vbar(
        x='plot_pos',
        width=100000000,
        top='box_max',
        bottom='box_min',
        color='black',
        fill_color=None,
        source=bokeh.models.ColumnDataSource({col: multi_assocs[col] for col in multi_assocs.columns})
    )

    results_plot.square(
        x='plot_pos',
        y='association_p_value',
        size='box_size',
        color='black',
        fill_color=None,
        source=bokeh.models.ColumnDataSource({col: multi_assocs[col] for col in multi_assocs.columns})
    )

    results_plot.quad(
        top = multi_assocs['max'] + 5,
        bottom = multi_assocs['min'] - 5,
        left = multi_assocs['plot_pos'] - 10000000,
        right = multi_assocs['plot_pos'] + 10000000,
        color = 'black',
        fill_color = None
    )
    '''

    results_plot.multi_line(
        [[pos, pos] for pos in df['plot_pos']],
        [[0, val] for val in df['association_p_value']],
        line_width=1,
        color=df['color']
    )
    circles = []
    for color, legend in (('blue', 'serum biomarker'), ('red', 'haematology phenotype')):
        color_df = df.filter(pl.col('color') == color)
        circles.append(results_plot.circle(
            'plot_pos',
            'association_p_value',
            size='size',
            fill_color='color',
            line_color='black',
            line_width=1,
            source=bokeh.models.ColumnDataSource({col: color_df[col] for col in color_df.columns}),
            legend_label = legend
        ))


    if args.ext == 'html':
        hover_tool = bokeh.models.tools.HoverTool(renderers = circles)
        hover_tool.tooltips = [
            ('phenotype', '@phenotype'),
            ('chrom', '@chr'),
            ('pos', '@pos'),
            ('repeat_unit', '@repeat_unit'),
            ('effect_direction', '@direction_of_association'),
            ('-log10(p_val)', '@association_p_value'),
            ('p_causal', '@pcausal'),
            ('relation_to_gene', '@relation_to_gene'),
            ('other_ethnic_effect_directions', '@other_ethnic_effect_directions')
        ]
        results_plot.add_tools(hover_tool)
        results_plot.toolbar.active_inspect = [hover_tool]

    '''
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    results_plot.add_layout(bokeh.models.Title(
        text=f"Plot creation date: {today}",
        align="right",
        text_font_size='18px'
    ))
    '''

    region_plots.export(results_plot, args.outfile, args.ext)

if __name__ == '__main__':
    main()
