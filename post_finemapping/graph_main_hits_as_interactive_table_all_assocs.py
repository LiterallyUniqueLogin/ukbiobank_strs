#!/usr/bin/env python3

import argparse
import os
import re

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
    parser.add_argument('e_splice_STRs_table')

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

    shaded_phenos = [
        'platelet_count',
        'platelet_crit',
        'mean_platelet_volume',
        'platelet_distribution_width',
        ###
        'cystatin_c',
        'creatinine',
        'urate',
        'urea',
        ###
        'igf_1',
        'shbg',
        'glucose',
        'calcium',
        'phosphate',
        ###
        'c_reactive_protein',
    ]
    pheno_blocks = {
        'red blood cell' : [
            'haematocrit',
            'haemoglobin_concentration',
            'red_blood_cell_count',
            'mean_corpuscular_volume',
            'mean_corpuscular_haemoglobin',
            'mean_sphered_cell_volume',
            'red_blood_cell_distribution_width',
            'glycated_haemoglobin',
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

    '''
    pheno_order = [
        ### red blood cell block
        'haematocrit',
        'haemoglobin_concentration',
        'red_blood_cell_count',
        'mean_corpuscular_volume',
        'mean_corpuscular_haemoglobin',
        'mean_sphered_cell_volume',
        'red_blood_cell_distribution_width',
        'glycated_haemoglobin',
        'mean_corpuscular_haemoglobin_concentration',
        ### platelet block
        'platelet_count',
        'platelet_crit',
        'mean_platelet_volume',
        'platelet_distribution_width',
        ### white blood cell block
        'eosinophil_count',
        'eosinophil_percent',
        'neutrophil_count',
        'neutrophil_percent',
        'lymphocyte_count',
        'lymphocyte_percent',
        'white_blood_cell_count',
        ### renal profile block
        'cystatin_c',
        'creatinine',
        'urate',
        'urea',
        ### liver profile block
        'gamma_glutamyltransferase',
        'alkaline_phosphatase',
        'total_bilirubin',
        'alanine_aminotransferase',
        'total_protein',
        'albumin',
        'aspartate_aminotransferase',
        ### endocrine profile block
        'igf_1',
        'shbg',
        'calcium',
        'glucose',
        'phosphate',
        ### lipid profile block
        'apolipoprotein_a',
        'hdl_cholesterol',
        'apolipoprotein_b',
        'ldl_cholesterol_direct',
        'cholesterol',
        'triglycerides',
        ### standalone
        'c_reactive_protein',
    ]
    '''

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

    for pair in [(23, 25), (25, 26), (30, 27), (32, 28), (33, 29), (29, 27), (37, 41), (38, 41), (39, 40), (42, 44), (47, 45), (58, 59), (61, 62), (50, 55), (51, 57), (52, 56), (50, 54), (51, 53), (70, 73), (84, 85),(85, 87), (87, 88)]:
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
    print(strs.dtype)
    _, idxs = np.unique(strs, return_index=True)
    unique_stable_strs = strs[np.sort(idxs)].flatten()
    unique_stable_strs = list(np.char.partition(np.array(unique_stable_strs, dtype=str), '_')[:, 2])
    # two genes each containing two unique hits
    unique_stable_strs[unique_stable_strs.index('CCDC26')] += ' #1'
    unique_stable_strs[unique_stable_strs.index('CCDC26')] += ' #2'
    unique_stable_strs[unique_stable_strs.index('TFDP2')] += ' #1'
    unique_stable_strs[unique_stable_strs.index('TFDP2')] += ' #2'

    print(unique_stable_strs)

    plot_width = 2500
    results_plot = bokeh.plotting.figure(
        width=plot_width,
        height=1400,
        x_axis_label='hits (containing gene, or position and nearest gene if intergenic)',
        y_axis_label='phenotypes',
        #y_range=bokeh.models.FactorRange(*[(group, pheno) for group in pheno_blocks for pheno in pheno_blocks[group]][::-1]),
        y_range=[pheno.replace('_', ' ') for pheno in pheno_order[::-1]],
        x_range=unique_stable_strs,#(0,94),
        tools='reset, save'
    )
    results_plot.xaxis.major_label_orientation = 1.5
    cds = bokeh.models.ColumnDataSource(dict(
        x=[len(x_coords)/2]*len(shaded_phenos), y=shaded_phenos, width=[len(x_coords)]*len(shaded_phenos), height=[1]*len(shaded_phenos), color=['grey']*len(shaded_phenos), alpha=['0.15']*len(shaded_phenos), line_color=[None]*len(shaded_phenos)
    ))
    results_plot.rect(
        x='x', y='y', width='width', height='height', color='color', alpha='alpha', line_color='line_color', source=cds
    )
    half_xs=np.arange(1, np.max(x_coords) + 1, 2)
    cds = bokeh.models.ColumnDataSource(dict(
        x=half_xs+0.5, y=[len(pheno_order)/2]*half_xs.shape[0], width=[1]*half_xs.shape[0], height=[len(pheno_order)+len(pheno_blocks)+20]*half_xs.shape[0], color=['grey']*half_xs.shape[0], alpha=['0.15']*half_xs.shape[0], line_color=[None]*half_xs.shape[0]
    ))
    results_plot.rect(
        x='x', y='y', width='width', height='height', color='color', alpha='alpha', line_color='line_color', source=cds
    )
    results_plot.xgrid.ticker = []

    str_colors = []
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

    def pheno_to_ycoords(phenos):
        return [pheno.replace('_', ' ') for pheno in phenos]
        #return [(group, pheno) for pheno in phenos for group in pheno_blocks if pheno in pheno_blocks[group]]

    undotted_loci = (df['finemapping'].to_numpy() == 'confidently') | (df['finemapping'].to_numpy() == 'not')

    undotted_triangle = (df['direction_of_association'].to_numpy() == '+') & undotted_loci
    results_plot.triangle(
        (x_coords + 0.5)[undotted_triangle],
        pheno_to_ycoords(df['phenotype'].to_numpy()[undotted_triangle]),
        #[(group, pheno) for pheno in df['phenotype'].to_numpy()[undotted_triangle] for group in pheno_blocks if pheno in pheno_blocks[group]],
        #df['phenotype'].to_numpy()[undotted_triangle],
        size=np.sqrt(df['p_val'].to_numpy())[undotted_triangle]*2,
        fill_alpha=df.select(pl.when(pl.col('finemapping') == 'confidently').then(1).otherwise(0)).to_numpy().flatten()[undotted_triangle],
        color='black'#str_colors[undotted_triangle],
    )
    undotted_inverted_triangle = (df['direction_of_association'].to_numpy() == '-') & undotted_loci
    results_plot.inverted_triangle(
        (x_coords + 0.5)[undotted_inverted_triangle],
        pheno_to_ycoords(df['phenotype'].to_numpy()[undotted_inverted_triangle]),
        #[(group, pheno) for pheno in df['phenotype'].to_numpy()[undotted_inverted_triangle] for group in pheno_blocks if pheno in pheno_blocks[group]],
        size=np.sqrt(df['p_val'].to_numpy())[undotted_inverted_triangle]*2,
        fill_alpha=df.select(pl.when(pl.col('finemapping') == 'confidently').then(1).otherwise(0)).to_numpy().flatten()[undotted_inverted_triangle],
        color='black'#str_colors[undotted_inverted_triangle],
    )

    dotted_triangle = (df['direction_of_association'].to_numpy() == '+') & ~undotted_loci
    results_plot.triangle_dot(
    #results_plot.triangle(
        (x_coords + 0.5)[dotted_triangle],
        pheno_to_ycoords(df['phenotype'].to_numpy()[dotted_triangle]),
        #[(group, pheno) for pheno in df['phenotype'].to_numpy()[dotted_triangle] for group in pheno_blocks if pheno in pheno_blocks[group]],
        #df['phenotype'].to_numpy()[dotted_triangle],
        size=np.sqrt(df['p_val'].to_numpy())[dotted_triangle]*2,
        fill_alpha=[0]*np.sum(dotted_triangle),
        color='black',#str_colors[dotted_triangle],
    )
    dotted_inverted_triangle = (df['direction_of_association'].to_numpy() == '-') & ~undotted_loci
    results_plot.triangle_dot(
    #results_plot.triangle(
        (x_coords + 0.5)[dotted_inverted_triangle],
        pheno_to_ycoords(df['phenotype'].to_numpy()[dotted_inverted_triangle]),
        #[(group, pheno) for pheno in df['phenotype'].to_numpy()[dotted_inverted_triangle] for group in pheno_blocks if pheno in pheno_blocks[group]],
        size=np.sqrt(df['p_val'].to_numpy())[dotted_inverted_triangle]*2,
        fill_alpha=[0]*np.sum(dotted_inverted_triangle),
        color='black',#str_colors[dotted_inverted_triangle],
        angle=np.pi
    )

    other_ethnicities = ['Black',  'South Asian', 'Chinese', 'Irish', 'White Other']
    results_topper = bokeh.plotting.figure(
        width=plot_width,
        height=300,
        y_range=bokeh.models.FactorRange(
            *[f'{ethnicity} replication' for ethnicity in other_ethnicities],
            'polyA',
            'polyAC',
            'polyCCG',
            'number of common alleles'
        ),
        x_range=results_plot.x_range,
        toolbar_location=None
    )
    results_topper.xgrid.ticker = []
    results_topper.xaxis.ticker = []
    cds = bokeh.models.ColumnDataSource(dict(
        x=half_xs+0.5, y=[7.5/2]*half_xs.shape[0], width=[1]*half_xs.shape[0], height=[7.5]*half_xs.shape[0], color=['grey']*half_xs.shape[0], alpha=['0.15']*half_xs.shape[0], line_color=[None]*half_xs.shape[0]
    ))
    results_topper.rect(
        x='x', y='y', width='width', height='height', color='color', alpha='alpha', line_color='line_color', source=cds
    )

    for ethnicity_num, ethnicity in enumerate(other_ethnicities):
        replicates = df.select((
            (pl.col('finemapping') == 'confidently') &
            (pl.col('other_ethnicity_effect_directions').str.split_exact(",", ethnicity_num+1).struct.field(f'field_{ethnicity_num}').str.strip() == pl.col('direction_of_association')) &
            (pl.col('other_ethnicity_association_p_values').str.split_exact(",", ethnicity_num+1).struct.field(f'field_{ethnicity_num}').str.strip().cast(float)*pl.col('n_confident_assocs') <= 0.05)
        ).alias('out'))['out'].to_numpy()
        results_topper.circle(
            (x_coords + 0.5)[replicates],
            [f'{ethnicity} replication']*np.sum(replicates),
            color='black',
        )
    for repeat_unit in 'A', 'AC', 'CCG':
        selection = df['repeat_unit'].to_numpy() == repeat_unit
        results_topper.circle(
            (x_coords + 0.5)[selection],
            [f'poly{repeat_unit}'] * np.sum(selection),
            color='black',
        )
    indices = []
    for index, coord in enumerate(x_coords):
        if index == list(x_coords).index(coord):
            indices.append(index)
    assert len(indices) == max(x_coords) + 1
    unique_stable_strs = list(np.char.partition(np.array(unique_stable_strs, dtype=str), '_')[:, 2])
    cds = bokeh.models.ColumnDataSource(dict(
        x=x_coords[indices]+0.5,
        y=['number of common alleles']*(int(np.max(x_coords)) + 1),
        width=[1]*(int(np.max(x_coords)) + 1),
        height=[1]*(int(np.max(x_coords)) + 1),
        color=['black']*(int(np.max(x_coords)) + 1),
        alpha=df['n_common_alleles'].to_numpy()[indices]/np.max(df['n_common_alleles'].to_numpy()),
        line_color=[None]*(int(np.max(x_coords)) + 1)
    ))
    results_topper.rect(
        x='x', y='y', width='width', height='height', color='color', alpha='alpha', line_color='line_color', source=cds
    )

    bokeh.io.export_png(bokeh.layouts.column(results_topper, results_plot), filename=args.outfile)

if __name__ == '__main__':
    main()
