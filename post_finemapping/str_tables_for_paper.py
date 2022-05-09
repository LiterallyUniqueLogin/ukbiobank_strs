#!/usr/bin/env python3

import ast
import json
import os

import pandas as pd
import polars as pl

import annotation_utils
import phenotypes

ukb = os.environ['UKB']

other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']

#just for testing
#phenotypes.phenotypes_in_use = ['platelet_distribution_width', 'platelet_count']
#phenotypes.phenotypes_in_use = ['albumin', 'calcium']#platelet_distribution_width', 'platelet_count']
#phenotypes.phenotypes_in_use = ['igf_1']

def dosages_to_frequencies(dosage_dict_str):
    dosages = ast.literal_eval(dosage_dict_str)
    # drop zero alleles
    dosages = { k: v for k, v in dosages.items() if v != 0 }
    total_dosage = sum(dosages.values())
    return json.dumps({ k: f'{v/total_dosage*100:.2f}%' for k,v in dosages.items() }).replace('"', '')

finemapping_dfs = []
for phenotype in phenotypes.phenotypes_in_use:
    df = pl.scan_csv(
        f'{ukb}/post_finemapping/intermediate_results/finemapping_all_concordance_{phenotype}.tab',
        sep='\t',
        dtypes={
            **{f'{ethnicity}_p_val': float for ethnicity in other_ethnicities},
            **{f'{ethnicity}_coeff': float for ethnicity in other_ethnicities},
            **{f'{ethnicity}_se': float for ethnicity in other_ethnicities}
        }
    ).filter('is_STR')
    fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    with open(fname) as tsv:
        header = tsv.readline().strip()
    assoc_df = pl.scan_csv(
        fname,
        sep='\t',
        skip_rows=1,
        has_header=False,
        with_column_names = lambda _: header.replace('0.05_significance_CI', 'foo', 1).replace('5e-8_significance_CI', 'bar', 1).split('\t') # these duplicate column names won't be used anyway
    ).select([
        'chrom',
        'pos',
        pl.col('subset_total_per_allele_dosages').alias('white_brit_allele_dosages')
    ])
    df = df.join(
        assoc_df,
        how='left',
        on=['chrom', 'pos']
    )
    for ethnicity in other_ethnicities:
        fname = f'{ukb}/association/results_finemapped_only/{ethnicity}/{phenotype}/my_str/results.tab'
        with open(fname) as tsv:
            header = tsv.readline().strip()
        assoc_df = pl.scan_csv(
            fname,
            sep='\t',
            skip_rows=1,
            has_header=False,
            with_column_names = lambda _: header.replace('0.05_significance_CI', 'foo', 1).replace('5e-8_significance_CI', 'bar', 1).split('\t') # these duplicate column names won't be used anyway
        ).select([
            'chrom',
            'pos',
            pl.col('subset_total_per_allele_dosages').alias(f'{ethnicity}_allele_dosages')
        ])
        df = df.join(
            assoc_df,
            how='left',
            on=['chrom', 'pos']
        )
    finemapping_dfs.append(df.collect())
finemapping_results = pl.concat(finemapping_dfs).rename({'pos': 'snpstr_pos'})

finemapping_results = finemapping_results.filter(
    (pl.col('p_val') <= 5e-8) &
    (
        ((pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 0)) | (pl.col('finemap_pip') >= 0.8)
    ).any().over(['chrom', 'snpstr_pos'])
)

pos_table = pl.read_csv(
    f'{ukb}/snpstr/flank_trimmed_vcf/vars.tab',
    sep='\t'
)

finemapping_results = finemapping_results.join(
    pos_table,
    how='left',
    on=['chrom', 'snpstr_pos']
)

repeat_units = pl.read_csv(
    f'{ukb}/snpstr/repeat_units.tab',
    sep='\t'
)

finemapping_results = finemapping_results.join(
    repeat_units,
    how='left',
    on=['chrom', 'snpstr_pos']
)

concordance_cols = pl.read_csv(
    f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_white_blood_cell_count.tab',
    sep='\t',
    n_rows=1
).columns
concordance_results = pl.concat([
    pl.read_csv(
        f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab',
        sep='\t',
        dtypes={col: (float if 'cs' not in col else int) for col in concordance_cols if 'finemap' in col or 'susie' in col or 'p_val' in col}
    ) for phenotype in phenotypes.phenotypes_in_use
    if not os.path.exists(f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab.empty')
]).filter('is_STR').with_column(pl.col('pos').alias('snpstr_pos'))

finemapping_results = finemapping_results.join(
    concordance_results,
    how='left',
    on=['phenotype', 'chrom', 'snpstr_pos']
).with_columns([
    pl.when(pl.col('susie_alpha').is_null()).then(None).when(pl.col('susie_cs') >= 0).then(pl.col('susie_alpha')).otherwise(0).alias('susie_CP'),
    pl.when(pl.col('susie_alpha_hardcall').is_null()).then(None).when(pl.col('susie_cs_hardcall') >= 0).then(pl.col('susie_alpha_hardcall')).otherwise(0).alias('susie_CP_best_guess_genotypes'),
    pl.when(pl.col('susie_alpha_ratio').is_null()).then(None).when(pl.col('susie_cs_ratio') >= 0).then(pl.col('susie_alpha_ratio')).otherwise(0).alias('susie_CP_prior_snps_over_strs'),
    pl.col('finemap_pip').alias('finemap_CP'),
    pl.col('finemap_pip_p_thresh').alias('finemap_CP_pval_thresh_5e-4'),
    pl.col('finemap_pip_mac').alias('finemap_CP_mac_thresh_100'),
    pl.col('finemap_pip_prior_std_derived').alias('finemap_CP_prior_effect_size_0.05%'),
    pl.col('finemap_pip_total_prob').alias('finemap_CP_prior_4_signals'),
    pl.col('finemap_pip_conv_tol').alias('finemap_CP_stopping_thresh_1e-4'),
    pl.col('finemap_pip_ratio').alias('finemap_CP_prior_snps_over_strs'),
    pl.col('finemap_pip_prior_std_low').alias('finemap_CP_prior_effect_size_0.0025%'),
    pl.col('pos').alias('start_pos')
])

# move to pandas for some code that's easier to write imperatively
finemapping_results = finemapping_results.to_pandas()
#print(finemapping_results.info(verbose=True))

# handle gene & transcript annotations
print("Loading gene and transcript annotations ...", flush=True)
#nearby_exon_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/nearby_exon_d_1000')
#closest_gene_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/closest_gene', distance=True)
#nearby_gene_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/nearby_gene_d_100000')

print("Loading high level intersections...", flush=True)
#transcript_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_transcript', bp_overlap=True)
gene_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_gene', bp_overlap=True)
exon_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_exon', bp_overlap=True)

print("Loading low level intersections...", flush=True)
CDS_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_CDS', bp_overlap=True)
five_prime_UTR_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_five_prime_UTR', bp_overlap=True)
three_prime_UTR_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_three_prime_UTR', bp_overlap=True)
UTR_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_UTR', bp_overlap=True)

sub_exon_types = pd.concat([CDS_intersect_merge, five_prime_UTR_intersect_merge, three_prime_UTR_intersect_merge, UTR_intersect_merge])

nrows = finemapping_results.shape[0]

#nearby_exons = ['none']*nrows
#closest_gene = [None]*nrows
#nearby_genes = ['none']*nrows
relation_to_gene = ['']*nrows
#transcribed = ['untranscribed']*nrows

print("Handling gene and transcript annotations ...", flush=True)
for idx in range(nrows):
    chrom = finemapping_results['chrom'][idx]
    start_pos = finemapping_results['start_pos'][idx]
    end_pos = finemapping_results['end_pos'][idx]

    '''
    print(sub_exon_types[['start_pos', 'STR_pos', 'STR_ID', 'annotation_feature_type', 'annotation_pos', 'annotation_end_pos']])
    print(sub_exon_types[['start_pos', 'STR_pos', 'STR_ID', 'annotation_feature_type', 'annotation_pos', 'annotation_end_pos']].dtypes)
    print(exon_intersect_merge[['start_pos', 'STR_pos', 'STR_ID', 'annotation_feature_type', 'annotation_pos', 'annotation_end_pos']])
    print(exon_intersect_merge[['start_pos', 'STR_pos', 'STR_ID', 'annotation_feature_type', 'annotation_pos', 'annotation_end_pos']].dtypes)
    for line in nearby_exon_merge[
        (nearby_exon_merge['chrom'] == chrom) & (nearby_exon_merge['STR_pos'] == start_pos)
    ].itertuples():
        if nearby_exons[idx] == 'none':
            nearby_exons[idx] = ''
        else:
            nearby_exons[idx] += ','
        nearby_exons[idx] += annotation_utils.distance(
            start_pos, end_pos, line.annotation_pos, line.annotation_end_pos
        ) + ":"
        nearby_exons[idx] += annotation_utils.get_relation(
            start_pos, end_pos, line.annotation_pos, line.annotation_end_pos, line.annotation_strand
        ) + ":"
        nearby_exons[idx] += str(line.annotation_pos) + ":"
        nearby_exons[idx] += str(line.annotation_end_pos) + ":"
        nearby_exons[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name') + ":"
        nearby_exons[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_type')

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

    for line in nearby_gene_merge[
        (nearby_gene_merge['chrom'] == chrom) & (nearby_gene_merge['STR_pos'] == start_pos)
    ].itertuples():
        if nearby_genes[idx] == 'none':
            nearby_genes[idx] = ''
        else:
            nearby_genes[idx] += ','
        nearby_genes[idx] += annotation_utils.distance(
            start_pos, end_pos, line.annotation_pos, line.annotation_end_pos
        ) + ":"
        nearby_genes[idx] += annotation_utils.get_relation(
            start_pos, end_pos, line.annotation_pos, line.annotation_end_pos, line.annotation_strand
        ) + ":"
        nearby_genes[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name') + ":"
        nearby_genes[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_type')
    '''

    # dict from gene to types of intersections
    gene_intersections = {}
    #partial_overlap = False
    for line in gene_intersect_merge[
        (gene_intersect_merge['chrom'].astype(int) == chrom) & (gene_intersect_merge['STR_pos'].astype(int) == start_pos)
    ].itertuples():
        gene_name = annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name')
        gene_type = annotation_utils.get_gff_kvp(line.annotation_info, 'gene_type')
        if gene_name in gene_intersections:
            # already described this gene, just hope the description is the same
            continue
        gene_intersections[gene_name] = [
            gene_type,
            {p: False for p in range(start_pos, end_pos+1)}, # bps covered by gene annotation
            {p: False for p in range(start_pos, end_pos+1)}, # bps covered by exon annotations
            set(), # types of sub-exon intersections
            {p: False for p in range(start_pos, end_pos+1)} # bps covered by sub-exon annotations
        ]
        for p in gene_intersections[gene_name][1]:
            if line.annotation_pos <= p <= line.annotation_end_pos:
                gene_intersections[gene_name][1][p] = True
            '''
            else:
                partial_overlap = True
            '''

    if len(gene_intersections) == 0:
        relation_to_gene[idx] = 'intergenic'

    '''
    if partial_overlap:
        relation_to_gene[idx] = 'intergenic;'
    '''

    if len(gene_intersections) > 1:
        relation_to_gene[idx] += 'multigene;'

    for line in exon_intersect_merge[
        (exon_intersect_merge['chrom'].astype(int) == chrom) & (exon_intersect_merge['STR_pos'].astype(int) == start_pos)
    ].itertuples():
        gene_name = annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name')
        if gene_name not in gene_intersections:
            raise ValueError(f"STR {chrom} {start_pos} gene name {gene_name} not in list of overlapping genes {set(gene_intersections.keys())}")
        for p in gene_intersections[gene_name][2]:
            if line.annotation_pos <= p <= line.annotation_end_pos:
                gene_intersections[gene_name][2][p] = True

    for line in sub_exon_types[
        (sub_exon_types['chrom'].astype(int) == chrom) & (sub_exon_types['STR_pos'].astype(int) == start_pos)
    ].itertuples():
        gene_name = annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name')
        if gene_name not in gene_intersections:
            raise ValueError(f"STR {chrom} {start_pos} gene name {gene_name} not in list of overlapping genes {set(gene_intersections.keys())}")
        if gene_intersections[gene_name][2] == False:
            raise ValueError(f"STR {chrom} {start_pos} intersects an exon type but isn't in an exon!")
        gene_intersections[gene_name][3].add(line.annotation_feature_type)
        for p in gene_intersections[gene_name][4]:
            if line.annotation_pos <= p <= line.annotation_end_pos:
                if not gene_intersections[gene_name][2][p]:
                    raise ValueError(f"STR {chrom} {start_pos} at position {p} has sub-exon type {line.annotation_feature_type} but is not exonic!")
                gene_intersections[gene_name][4][p] = True

    first = True
    #print(gene_intersections)
    for gene_name, gene_data in gene_intersections.items():
        if not first:
            relation_to_gene[idx] += ';'
        first = False
        types = gene_data[3]
        for p, is_gene in gene_data[1].items():
            if is_gene and not gene_data[2][p]:
                types.add('intron')
                break
        for p, is_gene in gene_data[2].items():
            if is_gene and not gene_data[4][p]:
                types.add('exon_other')
                break
        relation_to_gene[idx] += ','.join(types) + ":"
        relation_to_gene[idx] += gene_name + ":" + gene_data[0]

    '''
    for line in transcript_intersect_merge[
        (transcript_intersect_merge['chrom'] == chrom) & (transcript_intersect_merge['STR_pos'] == start_pos)
    ].itertuples():
        if transcribed[idx] == 'untranscribed':
            transcribed[idx] = 'transcribed;'
        else:
            transcribed[idx] += ','
        transcribed[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'transcript_name') + ":"
        transcribed[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'transcript_type') + ":"
        transcribed[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'transcript_support_level')
    '''

'''
finemapping_results['nearby_exons'] = nearby_exons
finemapping_results['closest_gene'] = closest_gene
finemapping_results['nearby_genes'] = nearby_genes
'''
finemapping_results['relation_to_gene'] = relation_to_gene
#finemapping_results['transcribed'] = transcribed

# done with gene annotations, move back to polars
finemapping_results = pl.DataFrame(finemapping_results)

finemapping_results = finemapping_results.select([
    'phenotype',
    'chrom',
    'start_pos',
    'end_pos',
    pl.col('region').alias('finemapping_region'),
    pl.col('unit').alias('repeat_unit'),
    pl.col('white_brit_allele_dosages').apply(dosages_to_frequencies).alias('white_brit_allele_frequencies'),
    pl.col('p_val').alias('association_p_value'),
    pl.when(pl.col('coeff') > 0).then('+').otherwise('-').alias('direction_of_association'),
    'relation_to_gene',
    pl.when(
        (pl.col('susie_CP') >= 0.8) &
        (pl.col('finemap_CP') >= 0.8) &
        (pl.col('susie_CP_best_guess_genotypes') >= 0.8) &
        (pl.col('finemap_CP_pval_thresh_5e-4') >= 0.8) &
        (pl.col('finemap_CP_mac_thresh_100') >= 0.8) &
        (pl.col('finemap_CP_prior_effect_size_0.05%') >= 0.8) &
        (pl.col('finemap_CP_prior_4_signals') >= 0.8) &
        (pl.col('finemap_CP_stopping_thresh_1e-4') >= 0.8)
    ).then(
        'confidently'
    ).when(
        (pl.col('susie_CP') >= 0.8) &
        (pl.col('finemap_CP') >= 0.8)
    ).then(
        'doubly'
    ).when(
        (pl.col('susie_CP') >= 0.8) |
        (pl.col('finemap_CP') >= 0.8)
    ).then(
        'singly'
    ).otherwise(
        'not'
    ).alias('finemapping'),
    'susie_CP',
    'finemap_CP',
    'susie_CP_best_guess_genotypes',
    'finemap_CP_pval_thresh_5e-4',
    'finemap_CP_mac_thresh_100',
    'finemap_CP_prior_effect_size_0.05%',
    'finemap_CP_prior_4_signals',
    'finemap_CP_stopping_thresh_1e-4',
    'susie_CP_prior_snps_over_strs',
    'finemap_CP_prior_snps_over_strs',
    'finemap_CP_prior_effect_size_0.0025%',
    pl.sum([
        pl.col(f'{ethnicity}_p_val').cast(str) + pl.lit(', ') for ethnicity in other_ethnicities
    ]).str.replace(', $', '').alias('other_ethnicity_association_p_values'),
    pl.sum([
        pl.when(pl.col(f'{ethnicity}_p_val') > .05).then('NA').when(pl.col('coeff') > 0).then('+').otherwise('-') + pl.lit(', ')
        for ethnicity in other_ethnicities
    ]).str.replace(', $', '').alias('other_ethnicity_effect_directions'),
    *[pl.col(f'{ethnicity}_allele_dosages').apply(dosages_to_frequencies).alias(f'{ethnicity}_allele_frequencies') for ethnicity in other_ethnicities],
])

finemapping_results.write_csv(f'{ukb}/post_finemapping/results/singly_finemapped_strs_for_paper.tab', sep='\t')
finemapping_results.sort(['chrom', 'start_pos']).write_csv(f'{ukb}/post_finemapping/results/singly_finemapped_strs_sorted.tab', sep='\t')

confident_results = finemapping_results.filter(
    (pl.col('finemapping') == 'confidently').any().over(['chrom', 'start_pos']) &
    (pl.col('association_p_value') <= 1e-10)
)
confident_results.write_csv(f'{ukb}/post_finemapping/results/confidently_finemapped_strs_for_paper.tab', sep='\t')
confident_results.sort(['chrom', 'start_pos']).write_csv(f'{ukb}/post_finemapping/results/confidently_finemapped_strs_sorted.tab', sep='\t')
