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
    ).filter(
        pl.col('finemapping') == 'confidently'
    )

    shaded_phenos = [
        'platelet_count',
        'platelet_crit',
        'mean_platelet_volume',
        'platelet_distribution_width',
        'urate',
        'creatinine',
        'ifg_1',
        'shbg',
        'apolipoprotein_a',
        'triglycerides',
    ]
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
        ### platelet block
        'platelet_count',
        'platelet_crit',
        'mean_platelet_volume',
        'platelet_distribution_width',
        ### white blood cell block
        'eosinophil_count',
        'eosinophil_percent',
        'lymphocyte_count',
        'lymphocyte_percent',
        'neutrophil_count',
        'white_blood_cell_count',
        ### renal profile block
        'urate',
        'creatinine',
        ### liver profile block
        'gamma_glutamyltransferase',
        'alkaline_phosphatase',
        'aspartate_aminotransferase',
        'total_bilirubin',
        ### encodcrine profile block
        'igf_1',
        'shbg',
        ### unknown
        'calcium',
        'c_reactive_protein',
        ### lipid profile block
        'apolipoprotein_a',
        'triglycerides',
        ### block
        'total_protein',
    ]

    pheno_indices = np.where(df['phenotype'].to_numpy()[:, None] == np.array(pheno_order)[None, :])[1]

    df = df.with_column(
        pl.Series(pheno_indices).alias('pheno_indices')
    ).with_row_count().with_column(
        (pl.col('pheno_indices') + pl.col('row_nr')/1e5).min().over(['chrom', 'start_pos']).alias('min_assoc_pheno_index')
    ).with_column(
        pl.when(
            pl.col('association_p_value') == 0
        ).then(300).otherwise(
            -pl.col('association_p_value').log10()
        ).alias('p_val')
    )

    mapi = df['min_assoc_pheno_index'].to_numpy().copy()
    for i in range(len(pheno_order)):
        if i not in np.floor(mapi):
            mapi[mapi > i] -= 1
    x_coords = np.floor(mapi).astype(int)
    for i in np.unique(mapi):
        num_dups = np.unique(mapi[(mapi < i) & (np.floor(i) == np.floor(mapi))]).shape[0]
        if num_dups > 0:
            x_coords[np.floor(mapi) > np.floor(i)] += 1
            x_coords[mapi == i] += num_dups
    '''
    for i in range(max(mapi) + 1):
        if i not in mapi:
            print('missing from init', i)
            exit()
    '''
    for i in range(max(x_coords) + 1):
        if i not in x_coords:
            print(i)
            exit()

    for pair in [(14, 15), (23, 26), (42,44), (50, 55), (51, 56), (52, 57), (60, 61)]:#[(6, 10)]:
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
    unique_stable_strs = list(np.char.partition(np.array(unique_stable_strs, dtype=str), '_')[:, 2])
    # two genes each containing two unique hits
    unique_stable_strs[unique_stable_strs.index('CCDC26')] += ' (1)'
    unique_stable_strs[unique_stable_strs.index('CCDC26')] += ' (2)'
    unique_stable_strs[unique_stable_strs.index('TFDP2')] += ' (1)'
    unique_stable_strs[unique_stable_strs.index('TFDP2')] += ' (2)'
    results_plot = bokeh.plotting.figure(
        width=2500,
        height=600,
        x_axis_label='hits',
        y_axis_label='phenotypes',
        y_range=pheno_order[::-1],
        x_range=unique_stable_strs,
        tools='reset, save'
    )
    results_plot.xaxis.major_label_orientation = 1.5
    cds = bokeh.models.ColumnDataSource(dict(
        x=[len(x_coords)/2]*len(shaded_phenos), y=shaded_phenos, width=[len(x_coords)]*len(shaded_phenos), height=[1]*len(shaded_phenos), color=['grey']*len(shaded_phenos), alpha=['0.15']*len(shaded_phenos), line_color=[None]*len(shaded_phenos)
    ))
    results_plot.rect(
        x='x', y='y', width='width', height='height', color='color', alpha='alpha', line_color='line_color', source=cds
    )
    strs, counts = np.unique(x_coords, return_counts = True)
    multi_pheno_strs = strs[counts > 1]
    results_plot.xgrid.ticker = multi_pheno_strs + 0.5
    results_plot.xgrid.grid_line_color = 'grey'
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
        color=str_colors[df['direction_of_association'].to_numpy() == '+']
    )
    results_plot.inverted_triangle(
        (x_coords + 0.5)[df['direction_of_association'].to_numpy() == '-'],
        df['phenotype'].to_numpy()[df['direction_of_association'].to_numpy() == '-'],
        size=np.sqrt(df['p_val'].to_numpy())[df['direction_of_association'].to_numpy() == '-']*2,
        color=str_colors[df['direction_of_association'].to_numpy() == '-']
    )

    bokeh.io.export_png(results_plot, filename=args.outfile)

if __name__ == '__main__':
    main()
