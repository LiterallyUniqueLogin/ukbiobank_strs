#!/usr/bin/env python3

import argparse
import re

import bokeh.io
import bokeh.models
import bokeh.plotting
import numpy as np
import polars as pl

import phenotypes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile')
    parser.add_argument('hits_table')
    parser.add_argument('e_splice_STRs_table')

    args = parser.parse_args()

    df = pl.read_csv(
        args.hits_table,
        sep='\t',
    )

    shaded_phenos = [
        'platelet_count',
        'platelet_crit',
        'mean_platelet_volume',
        'platelet_distribution_width',
        'neutrophil_count',
        'neutrophil_percent',
        'white_blood_cell_count',
        'apolipoprotein_a',
        'hdl_cholesterol',
        'cholesterol',
        'alkaline_phosphatase',
        'aspartate_aminotransferase',
        'c_reactive_protein',
        'shbg',
        'alanine_aminotransferase',
        'phosphate',
    ]
    pheno_order = [
        ### block
        'haematocrit',
        'haemoglobin_concentration',
        'red_blood_cell_count',
        'mean_corpuscular_volume',
        'mean_corpuscular_haemoglobin',
        'mean_sphered_cell_volume',
        'red_blood_cell_distribution_width',
        'glycated_haemoglobin',
        'mean_corpuscular_haemoglobin_concentration',
        ### block
        'platelet_count',
        'platelet_crit',
        'mean_platelet_volume',
        'platelet_distribution_width',
        ### block
        'eosinophil_count',
        'eosinophil_percent',
        'neutrophil_count',
        'neutrophil_percent',
        'lymphocyte_count',
        'lymphocyte_percent',
        'white_blood_cell_count',
        ### block
        'cystatin_c',
        'urate',
        'creatinine',
        'urea',
        ### block
        'apolipoprotein_a',
        'hdl_cholesterol',
        'apolipoprotein_b',
        'ldl_cholesterol_direct',
        'cholesterol',
        ### block
        'albumin',
        'total_protein',
        ### unknown
        'alkaline_phosphatase',
        'triglycerides',
        'aspartate_aminotransferase',
        'calcium',
        'c_reactive_protein',
        'igf_1',
        'shbg',
        'total_bilirubin',
        'alanine_aminotransferase',
        'glucose',
        'phosphate',
        'gamma_glutamyltransferase',
    ]

    pheno_indices = np.where(df['phenotype'].to_numpy()[:, None] == np.array(pheno_order)[None, :])[1]
    for phenotype in pheno_order:
        if not phenotype in phenotypes.phenotypes_in_use:
            print(phenotype)
            assert False

    df = df.with_column(
        pl.Series(pheno_indices).alias('pheno_indices')
    ).with_row_count().with_column(
        (pl.col('pheno_indices') + pl.col('row_nr')/1e5 + pl.when(pl.col('finemapping') != 'confidently').then(100000).otherwise(0)).min().over(['chrom', 'start_pos']).alias('min_assoc_pheno_index')
    ).with_column(
        pl.when(
            pl.col('association_p_value') == 0
        ).then(300).otherwise(
            -pl.col('association_p_value').log10()
        ).alias('p_val')
    )

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

    for pair in [(23, 25), (25, 26), (30, 27), (32, 28), (33, 29), (42, 44), (58, 59), (61, 62), (50, 55), (51, 57), (52, 56), (50, 54), (51, 53)]:
        x_coords[x_coords == pair[0]] = 1000000
        x_coords[x_coords == pair[1]] = pair[0]
        x_coords[x_coords == 1000000] = pair[1]

    strs = df[
        pl.when(
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
            pl.col('chrom').cast(str)+':'+pl.col('start_pos').cast(str)
        ).when(
            pl.col('relation_to_gene').str.contains('multigene')
        ).then(
            pl.col('relation_to_gene').str.extract(r"\w[^{:]*:([^:]*):", 1) + ' & ' +
            pl.col('relation_to_gene').str.extract(r"\w[^{:;]*;[^:]*:([^:]*):", 1)
        ).otherwise(
            pl.col('relation_to_gene').str.extract(r"\w[^{:]*:([^:]*):", 1)
        )
    ].to_numpy()[np.argsort(x_coords)]
    _, idxs = np.unique(strs, return_index=True)
    unique_stable_strs = strs[np.sort(idxs)].flatten()
    results_plot = bokeh.plotting.figure(
        width=2500,
        height=900,
        x_axis_label='hits',
        y_axis_label='phenotypes',
        y_range=pheno_order[::-1],
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
        x=half_xs+0.5, y=[len(pheno_order)/2]*half_xs.shape[0], width=[1]*half_xs.shape[0], height=[len(pheno_order)]*half_xs.shape[0], color=['grey']*half_xs.shape[0], alpha=['0.15']*half_xs.shape[0], line_color=[None]*half_xs.shape[0]
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

    results_plot.triangle(
        (x_coords + 0.5)[df['direction_of_association'].to_numpy() == '+'],
        df['phenotype'].to_numpy()[df['direction_of_association'].to_numpy() == '+'],
        size=np.sqrt(df['p_val'].to_numpy())[df['direction_of_association'].to_numpy() == '+']*2,
        fill_alpha=df.select(pl.when(pl.col('finemapping') == 'confidently').then(1).otherwise(0)).to_numpy().flatten()[df['direction_of_association'].to_numpy() == '+'],
        color=str_colors[df['direction_of_association'].to_numpy() == '+'],
    )
    results_plot.inverted_triangle(
        (x_coords + 0.5)[df['direction_of_association'].to_numpy() == '-'],
        df['phenotype'].to_numpy()[df['direction_of_association'].to_numpy() == '-'],
        size=np.sqrt(df['p_val'].to_numpy())[df['direction_of_association'].to_numpy() == '-']*2,
        fill_alpha=df.select(pl.when(pl.col('finemapping') == 'confidently').then(1).otherwise(0)).to_numpy().flatten()[df['direction_of_association'].to_numpy() == '-'],
        color=str_colors[df['direction_of_association'].to_numpy() == '-'],
    )

    bokeh.io.export_png(results_plot, filename=args.outfile)

if __name__ == '__main__':
    main()
