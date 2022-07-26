#!/usr/bin/env python3

import argparse
import os

import bokeh.io
import bokeh.models
import bokeh.plotting
import bokeh.layouts
import numpy as np
import polars as pl

import annotation_utils
import phenotypes

ukb = os.environ['UKB']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile')
    parser.add_argument('hits_table')
    parser.add_argument('qtl_STRs_table')

    args = parser.parse_args()

    df = pl.read_csv(
        args.hits_table,
        sep='\t',
    ).filter(pl.col('association_p_value') <= 1e-10)

    closest_gene_merge = annotation_utils.get_merged_annotations(
        df.with_column(pl.col('start_pos').alias('pos')).to_pandas(),
        f'{ukb}/side_analyses/str_annotations/closest_gene',
        distance=True
    )
    closest_gene = [None]*df.shape[0]

    nrows = df.shape[0]
    for idx in range(nrows):
        chrom = df['chrom'][idx]
        start_pos = df['start_pos'][idx]
        end_pos = df['end_pos'][idx]
        for line in closest_gene_merge[
            (closest_gene_merge['chrom'] == chrom) & (closest_gene_merge['STR_pos'] == start_pos)
        ].itertuples():
            if closest_gene[idx] is not None:
                closest_gene[idx] += ','
            else:
                closest_gene[idx] = ''
            closest_gene[idx] += str(line.annotation_distance) + ":"
            closest_gene[idx] += annotation_utils.get_relation(
                start_pos, end_pos, line.annotation_pos, line.annotation_end_pos, line.annotation_strand
            ) + ":"

            closest_gene[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name') + ":"
            closest_gene[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_type')

    df = df.with_column(pl.Series(closest_gene).alias('closest_gene'))

    pheno_blocks = {
        'red blood cell' : [
            'haematocrit',
            'haemoglobin_concentration',
            'red_blood_cell_count',
            'mean_corpuscular_volume',
            'mean_corpuscular_haemoglobin',
            'mean_sphered_cell_volume',
            'red_blood_cell_distribution_width',
            'mean_corpuscular_haemoglobin_concentration',
        ],
        'platelet' : [
            'platelet_count',
            'platelet_crit',
            'mean_platelet_volume',
            'platelet_distribution_width',
        ],
        'white blood cell' : [
            'eosinophil_count',
            'eosinophil_percent',
            'neutrophil_count',
            'neutrophil_percent',
            'lymphocyte_count',
            'lymphocyte_percent',
            'white_blood_cell_count',
        ],
        'renal' : [
            'cystatin_c',
            'creatinine',
            'urate',
            'urea',
        ],
        'liver' : [
            'gamma_glutamyltransferase',
            'alkaline_phosphatase',
            'total_bilirubin',
            'alanine_aminotransferase',
            'total_protein',
            'albumin',
            'aspartate_aminotransferase',
        ],
        'endocrine' : [
            'glycated_haemoglobin',
            'igf_1',
            'shbg',
            'calcium',
            'glucose',
            'phosphate',
        ],
        'lipid' : [
            'apolipoprotein_a',
            'hdl_cholesterol',
            'apolipoprotein_b',
            'ldl_cholesterol_direct',
            'cholesterol',
            'triglycerides',
        ],
        '' : ['c_reactive_protein']
    }
    pheno_order = sum(pheno_blocks.values(), [])
    shaded_phenos = sum([pheno_blocks[key] for key in list(pheno_blocks.keys())[1::2]], [])

    pheno_indices = np.where(df['phenotype'].to_numpy()[:, None] == np.array(pheno_order)[None, :])[1]
    for phenotype in pheno_order:
        if not phenotype in phenotypes.phenotypes_in_use:
            print(phenotype)
            assert False

    df = df.with_column(
        pl.Series(pheno_indices).alias('pheno_indices')
    ).with_row_count().with_column(
        (pl.col('pheno_indices') + pl.col('row_nr')/1e5 + pl.when(pl.col('finemapping') != 'confidently').then(100000).otherwise(0)).min().over(['chrom', 'start_pos']).alias('min_assoc_pheno_index')
    ).with_columns([
        pl.when(
            pl.col('association_p_value') == 0
        ).then(300).otherwise(
            -pl.col('association_p_value').log10()
        ).alias('p_val'),
        (pl.col('finemapping') == 'confidently').sum().over(['chrom', 'start_pos']).alias('n_confident_assocs'),
        pl.col('white_brit_allele_frequencies').apply(
            lambda dosage_dict_str: sum(float(part.split(' ')[-1]) >= 1 for part in dosage_dict_str.split('%')[:-1])
        ).alias('n_common_alleles')
    ])

    mapi = df['min_assoc_pheno_index'].to_numpy().copy()
    old_mapi = mapi.copy()
    for i in range(len(pheno_order)):
        if i not in np.floor(old_mapi):
            mapi[old_mapi > i] -= 1
    x_coords = np.floor(mapi).astype(int)
    for i in np.unique(mapi):
        num_dups = np.unique(mapi[(mapi < i) & (np.floor(i) == np.floor(mapi))]).shape[0]
        if num_dups > 0:
            x_coords[np.floor(mapi) > np.floor(i)] += 1
            x_coords[mapi == i] += num_dups
    for i in range(max(x_coords) + 1):
        if i not in x_coords:
            print(i)
            exit()

    for pair in [(21, 23), (23, 24), (25, 28), (26, 30), (27, 31), (25, 27), (35, 39), (36, 39), (37, 38), (40, 42), (40, 41), (43, 45), (56, 57), (59, 60), (48, 53), (49, 55), (50, 54), (49, 51), (68, 71), (84, 85),(85, 87), (87, 88)]:
        x_coords[x_coords == pair[0]] = 1000000
        x_coords[x_coords == pair[1]] = pair[0]
        x_coords[x_coords == 1000000] = pair[1]

    strs = df[
        pl.col('chrom').cast(str) + ':' + pl.col('start_pos').cast(str) + '_' + pl.when(
            pl.col('relation_to_gene').str.contains('protein_coding.*protein_coding')
        ).then(
            pl.col('relation_to_gene').str.extract("([^:]*):protein_coding.*:([^:]*):protein_coding", 1) + ' & ' +
            pl.col('relation_to_gene').str.extract("([^:]*):protein_coding.*:([^:]*):protein_coding", 2)
        ).when(
            pl.col('relation_to_gene').str.contains('protein_coding')
        ).then(
            pl.col('relation_to_gene').str.extract("([^:]*):protein_coding", 1)
        ).when(
            pl.col('relation_to_gene').str.contains('intergenic')
        ).then(
            pl.col('chrom').cast(str)+':'+pl.col('start_pos').cast(str)+' (' + pl.col('closest_gene').str.split_exact(',', 1).struct.field('field_0').str.split_exact(':', 2).struct.field('field_2') + ')'
        ).when(
            pl.col('relation_to_gene').str.contains('multigene')
        ).then(
            pl.col('relation_to_gene').str.extract(r"multigene;[^:]*:([^:]*):", 1) + '/' +
            pl.col('relation_to_gene').str.extract(r"multigene;[^;]*;[^:]*:([^:]*):", 1)
        ).otherwise(
            pl.col('relation_to_gene').str.extract(r"\w[^{:]*:([^:]*):", 1)
        )
    ].to_numpy()[np.argsort(x_coords)]
    _, idxs = np.unique(strs, return_index=True)
    unique_stable_strs = strs[np.sort(idxs)].flatten()
    unique_stable_strs = list(np.char.partition(np.array(unique_stable_strs, dtype=str), '_')[:, 2])
    # two genes each containing two unique hits
    unique_stable_strs[unique_stable_strs.index('CCDC26')] += ' #1'
    unique_stable_strs[unique_stable_strs.index('CCDC26')] += ' #2'
    unique_stable_strs[unique_stable_strs.index('TFDP2')] += ' #1'
    unique_stable_strs[unique_stable_strs.index('TFDP2')] += ' #2'

    plot_width = 2500
    results_plot = bokeh.plotting.figure(
        width=plot_width,
        height=1400,
        x_axis_label='STRs (containing gene, or position and nearest gene if intergenic)',
        y_axis_label='phenotypes',
        #y_range=bokeh.models.FactorRange(*[(group, pheno) for group in pheno_blocks for pheno in pheno_blocks[group]][::-1]),
        y_range=[pheno.replace('_', ' ') for pheno in pheno_order[::-1]],
        x_range=unique_stable_strs,#(0,94),
        toolbar_location=None,
        outline_line_color='black',
        output_backend='svg'
    )
    results_plot.axis.axis_label_text_font_size = '30px'
    results_plot.xaxis.major_label_orientation = 1.3
    cds = bokeh.models.ColumnDataSource(dict(
        x=[len(x_coords)/2]*len(shaded_phenos), y=[pheno.replace('_', ' ') for pheno in shaded_phenos], width=[len(x_coords)]*len(shaded_phenos), height=[1]*len(shaded_phenos), color=['grey']*len(shaded_phenos), alpha=['0.15']*len(shaded_phenos), line_color=[None]*len(shaded_phenos)
    ))
    results_plot.rect(
        x='x', y='y', width='width', height='height', color='color', alpha='alpha', line_color='line_color', source=cds
    )
    rect_steps = np.arange(1, np.max(x_coords) + 1, 2)
    cds = bokeh.models.ColumnDataSource(dict(
        x=rect_steps+0.5,
        y=[len(pheno_order)/2]*rect_steps.shape[0],
        width=[1]*rect_steps.shape[0],
        height=[len(pheno_order)+len(pheno_blocks)+20]*rect_steps.shape[0],
        color=['grey']*rect_steps.shape[0],
        alpha=['0.075']*rect_steps.shape[0],
        line_color=[None]*rect_steps.shape[0]
    ))
    results_plot.rect(
        x='x', y='y', width='width', height='height', color='color', alpha='alpha', line_color='line_color', source=cds
    )
    results_plot.xgrid.ticker = []

    '''
    with open(args.e_splice_STRs_table) as table:
        e_splice_lines = table.readlines()

    for loc, gene_rels in zip(
        df[pl.col('chrom').cast(str)+'_'+pl.col('start_pos').cast(str)].to_numpy().flatten(),
        df['relation_to_gene'].to_numpy().flatten()
    ):
        gene_rel_list = gene_rels.split(';')
        color = 'blue'
        for gene_rel in gene_rel_list:
            if gene_rel == 'intergenic':
                break
            if gene_rel == 'multigene':
                continue
            for line in e_splice_lines:
                if re.search(loc + '.*' + gene_rel.split(':')[1], line):
                    color = 'red'

        str_colors.append(color)
    str_colors = np.array(str_colors)
    '''

    def pheno_to_ycoords(phenos):
        return [pheno.replace('_', ' ') for pheno in phenos]
        #return [(group, pheno) for pheno in phenos for group in pheno_blocks if pheno in pheno_blocks[group]]

    fill_colors = np.array(['           ']*len(x_coords))
    fill_colors[:] = '#A13D34' #'#C52F22' #'grey'
    fill_colors[df['finemapping'].to_numpy() == 'confidently'] = '#401E1D' #'black'
    fill_colors[df['finemapping'].to_numpy() == 'not'] = '#BF9E9A' #'white'

    up_triangle = df['direction_of_association'].to_numpy() == '+'
    results_plot.triangle(
        (x_coords + 0.5)[up_triangle],
        pheno_to_ycoords(df['phenotype'].to_numpy()[up_triangle]),
        size=np.sqrt(df['p_val'].to_numpy())[up_triangle]*2,
        color=fill_colors[up_triangle],
        line_color=None
    )
    results_plot.inverted_triangle(
        (x_coords + 0.5)[~up_triangle],
        pheno_to_ycoords(df['phenotype'].to_numpy()[~up_triangle]),
        size=np.sqrt(df['p_val'].to_numpy())[~up_triangle]*2,
        fill_color=fill_colors[~up_triangle],
        line_color=None
    )

    other_ethnicities = ['Black',  'South Asian', 'Chinese', 'Irish', 'White Other']

    topper_circle_size = 8
    def get_topper(height, factors, color):
        topper = bokeh.plotting.figure(
            width=plot_width,
            height=height,
            y_range=bokeh.models.FactorRange(*factors),
            x_range=results_plot.x_range,
            toolbar_location=None,
            outline_line_color='black',
            output_backend='svg'
        )
        topper.xgrid.ticker = []
        topper.xaxis.ticker = []
        if color:
            cds = bokeh.models.ColumnDataSource(dict(
                x=rect_steps+0.5,
                y=[7.5/2]*rect_steps.shape[0],
                width=[1]*rect_steps.shape[0],
                height=[7.5]*rect_steps.shape[0],
                color=['grey']*rect_steps.shape[0],
                alpha=['0.075']*rect_steps.shape[0],
                line_color=[None]*rect_steps.shape[0]
            ))
            topper.rect(
                x='x', y='y', width='width', height='height', color='color', alpha='alpha', line_color='line_color', source=cds
            )
        return topper
    
    qtl_topper = get_topper(
        30, ['expression QTL'], True#, 'splice or isoform QTL'], True
    )

    qtl_STRs = pl.read_csv(args.qtl_STRs_table, sep='\t')

    eqtl_STR_locs = qtl_STRs.filter(~pl.col('p_vals_expression').is_null())['chrom_pos']
    eqtl_STRs = df[
        ('chr' + pl.col('chrom').cast(str) + '_' + pl.col('start_pos').cast(str)).is_in(eqtl_STR_locs)
    ].to_numpy().flatten()
    qtl_topper.circle(
        (x_coords + 0.5)[eqtl_STRs],
        ['expression QTL']*np.sum(eqtl_STRs),
        color='black',
        size=topper_circle_size
    )

    '''
    splice_iso_STR_locs = qtl_STRs.filter(~pl.col('p_vals_splice').is_null() | ~pl.col('p_vals_isoform').is_null())['chrom_pos']
    splice_iso_STRs = df[
        ('chr' + pl.col('chrom').cast(str) + '_' + pl.col('start_pos').cast(str)).is_in(splice_iso_STR_locs)
    ].to_numpy().flatten()
    qtl_topper.circle(
        (x_coords + 0.5)[splice_iso_STRs],
        ['splice or isoform QTL']*np.sum(splice_iso_STRs),
        color='black',
        size=topper_circle_size
    )
    '''


    replication_topper = get_topper(
        150, [f'{ethnicity} replication'.replace('Other', 'other') for ethnicity in other_ethnicities], True
    )

    for ethnicity_num, ethnicity in enumerate(other_ethnicities):
        replicates = df.select((
            (pl.col('finemapping') == 'confidently') &
            (pl.col('other_ethnicity_effect_directions').str.split_exact(",", ethnicity_num+1).struct.field(f'field_{ethnicity_num}').str.strip() == pl.col('direction_of_association')) &
            (pl.col('other_ethnicity_association_p_values').str.split_exact(",", ethnicity_num+1).struct.field(f'field_{ethnicity_num}').str.strip().cast(float)*pl.col('n_confident_assocs') <= 0.05)
        ).alias('out'))['out'].to_numpy()
        replication_topper.circle(
            (x_coords + 0.5)[replicates],
            [f'{ethnicity} replication'.replace('Other', 'other')]*np.sum(replicates),
            color='black',
            size=topper_circle_size
        )

    repeat_unit_topper = get_topper(90, ['polyA', 'polyAC', 'polyCCG'], True)
    for repeat_unit in 'A', 'AC', 'CCG':
        selection = df['repeat_unit'].to_numpy() == repeat_unit
        repeat_unit_topper.circle(
            (x_coords + 0.5)[selection],
            [f'poly{repeat_unit}'] * np.sum(selection),
            color='black',
            size=topper_circle_size
        )


    indices = []
    for index, coord in enumerate(x_coords):
        if index == list(x_coords).index(coord):
            indices.append(index)
    assert len(indices) == max(x_coords) + 1
    unique_stable_strs = list(np.char.partition(np.array(unique_stable_strs, dtype=str), '_')[:, 2])

    multiallelic_topper = get_topper(30, ['number of common alleles'], False)

    max_common_alleles = np.max(df['n_common_alleles'].to_numpy())
    cds = bokeh.models.ColumnDataSource(dict(
        x=x_coords[indices]+0.5,
        y=['number of common alleles']*(int(np.max(x_coords)) + 1),
        width=[1]*(int(np.max(x_coords)) + 1),
        height=[1]*(int(np.max(x_coords)) + 1),
        color=['black']*(int(np.max(x_coords)) + 1),
        alpha=df['n_common_alleles'].to_numpy()[indices]/max_common_alleles,
        line_color=[None]*(int(np.max(x_coords)) + 1)
    ))
    multiallelic_topper.rect(
        x='x', y='y', width='width', height='height', color='color', alpha='alpha', line_color='line_color', source=cds
    )

    multiallelic_scale = bokeh.plotting.figure(
        width=120,
        height=40*max_common_alleles,
        y_range=[0.5, max_common_alleles+0.5],
        y_axis_label='# common alleles',
        toolbar_location=None,
        outline_line_color='black',
        output_backend='svg'
    )
    multiallelic_scale.axis.axis_label_text_font_size = '26px'
    multiallelic_scale.xaxis.ticker = []
    multiallelic_scale.xgrid.ticker = []
    multiallelic_scale.yaxis.ticker = np.arange(1, max_common_alleles + 1)

    cds = bokeh.models.ColumnDataSource(dict(
        x=[0]*max_common_alleles,
        y=np.arange(1, max_common_alleles+1),
        width=[1]*max_common_alleles,
        height=[1]*max_common_alleles,
        color=['black']*max_common_alleles,
        alpha=np.arange(1, max_common_alleles+1)/max_common_alleles,
        line_color=[None]*max_common_alleles
    ))
    multiallelic_scale.rect(
        x='x', y='y', width='width', height='height',
        color='color', alpha='alpha', line_color='line_color',
        source=cds
    )

    scale_ps = [10, 20, 40, 80, 160, 300]
    str_scale_ps = [str(p) for p in scale_ps]
    p_val_scale = bokeh.plotting.figure(
        width=120,
        height=30*len(scale_ps) + 90,
        y_range=bokeh.models.FactorRange(*str_scale_ps),
        y_axis_label='-log10 p-value',
        toolbar_location=None,
        outline_line_color='black',
        output_backend='svg'
    )
    p_val_scale.axis.axis_label_text_font_size = '26px'
    p_val_scale.xaxis.ticker = []
    p_val_scale.grid.ticker = []
    p_val_scale.triangle(
        [0]*len(scale_ps),
        str_scale_ps,
        size=np.sqrt(scale_ps)*2,
        color='black'
    )

    total_fig = bokeh.layouts.row(bokeh.layouts.column(
            multiallelic_topper, repeat_unit_topper, replication_topper, qtl_topper, results_plot, 
        ),
        bokeh.layouts.column(multiallelic_scale, p_val_scale)
    )
    bokeh.io.export_png(
        total_fig,
        filename=args.outfile
    )

    bokeh.io.export_svg(
        total_fig,
        filename=args.outfile.replace('png', 'svg')
    )

if __name__ == '__main__':
    main()
