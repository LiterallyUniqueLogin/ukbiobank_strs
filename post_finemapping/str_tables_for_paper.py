#!/usr/bin/env python3

import argparse
import ast
import json

import pandas as pd
import polars as pl

import annotation_utils

other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']

parser = argparse.ArgumentParser()
parser.add_argument('--first-pass-finemapping-dfs', nargs='+')
parser.add_argument('--followup-finemapping-dfs', nargs='+')
parser.add_argument('--assoc-phenotypes', nargs='+')
parser.add_argument('--assocs', nargs='+')
parser.add_argument('--black-assocs', nargs='+')
parser.add_argument('--south-asian-assocs', nargs='+')
parser.add_argument('--chinese-assocs', nargs='+')
parser.add_argument('--irish-assocs', nargs='+')
parser.add_argument('--white-other-assocs', nargs='+')
parser.add_argument('--str-pos-table')
parser.add_argument('--str-pos-table-2')
parser.add_argument('--str-hg38-pos-table')
parser.add_argument('--str-t2t-pos-table')
parser.add_argument('--repeat-units-table')
parser.add_argument('--intersects-gene-annotation', nargs=22)
parser.add_argument('--intersects-exon-annotation', nargs=22)
parser.add_argument('--intersects-CDS-annotation', nargs=22)
parser.add_argument('--intersects-five-prime-UTR-annotation', nargs=22)
parser.add_argument('--intersects-three-prime-UTR-annotation', nargs=22)
parser.add_argument('--intersects-UTR-annotation', nargs=22)
parser.add_argument('--outdir')

args = parser.parse_args()

assert len(args.assoc_phenotypes) == \
        len(args.assocs) == \
        len(args.black_assocs) == \
        len(args.south_asian_assocs) == \
        len(args.chinese_assocs) == \
        len(args.irish_assocs) == \
        len(args.white_other_assocs)

assocs = dict(zip(args.assoc_phenotypes, args.assocs))
other_ethnicity_assocs = {
    'black': dict(zip(args.assoc_phenotypes, args.black_assocs)),
    'south_asian': dict(zip(args.assoc_phenotypes, args.south_asian_assocs)),
    'chinese': dict(zip(args.assoc_phenotypes, args.chinese_assocs)),
    'irish': dict(zip(args.assoc_phenotypes, args.irish_assocs)),
    'white_other': dict(zip(args.assoc_phenotypes, args.white_other_assocs))
}

def dosages_to_frequencies(dosage_dict_str):
    dosages = ast.literal_eval(dosage_dict_str)
    # drop zero alleles
    dosages = { k: v for k, v in dosages.items() if v != 0 }
    total_dosage = sum(dosages.values())
    return json.dumps({ k: f'{v/total_dosage*100:.2f}%' for k,v in dosages.items() }).replace('"', '')

finemapping_dfs = []
# f'{ukb}/post_finemapping/intermediate_results/finemapping_all_concordance_{phenotype}.tab',
for count, fname in enumerate(args.first_pass_finemapping_dfs):
    print(f'Loading finemapping df {count+1}/{len(args.first_pass_finemapping_dfs)}', flush=True)
    df = pl.DataFrame(pd.read_csv(
        fname,
        sep='\t',
        dtype={
            **{f'{ethnicity}_p_val': float for ethnicity in other_ethnicities},
            **{f'{ethnicity}_coeff': float for ethnicity in other_ethnicities},
            **{f'{ethnicity}_se': float for ethnicity in other_ethnicities}
        }
    )).filter('is_STR')
    phenotype = df['phenotype'][0] # all will be the same
    # need to add phenotype col
    # fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    assoc_df = pl.scan_csv(
        assocs[phenotype],
        sep='\t',
    ).select([
        'chrom',
        'pos',
        pl.col('subset_total_per_allele_dosages').alias('white_brit_allele_dosages')
    ])
    lazy_df = df.lazy().join(
        assoc_df,
        how='left',
        on=['chrom', 'pos']
    )
    for ethnicity in other_ethnicities:
        # need to add phenotype col
        # fname = f'{ukb}/association/results_finemapped_only/{ethnicity}/{phenotype}/my_str/results.tab'
        assoc_df = pl.scan_csv(
            other_ethnicity_assocs[ethnicity][phenotype],
            sep='\t',
        ).select([
            'chrom',
            'pos',
            pl.col('subset_total_per_allele_dosages').alias(f'{ethnicity}_allele_dosages')
        ])
        lazy_df = lazy_df.join(
            assoc_df,
            how='left',
            on=['chrom', 'pos']
        )
    finemapping_dfs.append(lazy_df.collect())
finemapping_results = pl.concat(finemapping_dfs).rename({'pos': 'snpstr_pos'})

finemapping_results = finemapping_results.filter(
    (pl.col('p_val') <= 5e-8) &
    (
        ((pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 0)) | (pl.col('finemap_pip') >= 0.8)
    ).any().over(['chrom', 'snpstr_pos'])
)

pos_table = pl.read_csv(
    args.str_pos_table, sep='\t'
).groupby(
    ['chrom', 'pos', 'end_pos']
).agg([pl.col('snpstr_pos').first()]) # drop duplicate snpstr poses

pos_table_2 = pl.read_csv(
    args.str_pos_table_2,
    sep='\t',
    has_header=False,
    new_columns=['chrom', 'pos', 'end_pos', 'ID']
).with_columns([
    pl.col('pos') + 1, # transform from bed to vcf coords
    pl.col('chrom').str.replace('chr', '').cast(int),
])
hg38_pos_table = pl.read_csv(
    args.str_hg38_pos_table,
    sep='\t',
    has_header=False,
    new_columns=['chrom', 'pos', 'end_pos', 'ID', 'drop']
).with_columns([
    pl.col('pos') + 1, # transform from bed to vcf coords
]).drop(['drop', 'chrom'])
t2t_pos_table = pl.read_csv(
    args.str_t2t_pos_table,
    sep='\t',
    has_header=False,
    new_columns=['chrom', 'pos', 'end_pos', 'ID', 'drop']
).with_columns([
    pl.col('pos') + 1, # transform from bed to vcf coords
]).drop(['drop', 'chrom'])

pos_table_2 = pos_table_2.join(
    hg38_pos_table,
    how = 'left',
    on=['ID'],
    suffix='_hg38',
).join(
    t2t_pos_table,
    how = 'left',
    on=['ID'],
    suffix='_t2t'
)

pos_table = pos_table_2.join(
    pos_table,
    on=['chrom', 'pos', 'end_pos']
).drop(['ID'])

finemapping_results = finemapping_results.join(
    pos_table,
    how='left',
    on=['chrom', 'snpstr_pos']
)

repeat_units = pl.read_csv(
    #f'{ukb}/snpstr/repeat_units.tab',
    args.repeat_units_table,
    sep='\t'
)

finemapping_results = finemapping_results.join(
    repeat_units,
    how='left',
    on=['chrom', 'snpstr_pos']
)

concordance_cols = pl.read_csv(
    #f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_white_blood_cell_count.tab',
    args.followup_finemapping_dfs[0],
    sep='\t',
    n_rows=1
).columns
concordance_results = pl.concat([
    pl.read_csv(
        #f'{ukb}/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_{phenotype}.tab',
        followup_finemapping_df,
        sep='\t',
        dtypes={col: (float if 'cs' not in col else int) for col in concordance_cols if 'finemap' in col or 'susie' in col or 'p_val' in col}
    ) for followup_finemapping_df in args.followup_finemapping_dfs
]).filter('is_STR').with_column(pl.col('pos').alias('snpstr_pos'))

finemapping_results = finemapping_results.join(
    concordance_results,
    how='left',
    on=['phenotype', 'chrom', 'snpstr_pos']
).with_columns([
    pl.when(pl.col('susie_alpha').is_null()).then(None).when(pl.col('susie_cs') >= 0).then(pl.col('susie_alpha')).otherwise(0).alias('susie_CP'),
    pl.when(pl.col('susie_alpha_best_guess').is_null()).then(None).when(pl.col('susie_cs_best_guess') >= 0).then(pl.col('susie_alpha_best_guess')).otherwise(0).alias('susie_CP_best_guess_genotypes'),
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

# handle gene & transcript annotations
print("Loading gene and transcript annotations ...", flush=True)

print("Loading high level intersections...", flush=True)
#gene_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_gene', bp_overlap=True)
#exon_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_exon', bp_overlap=True)
gene_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, args.intersects_gene_annotation, bp_overlap=True)
exon_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, args.intersects_exon_annotation, bp_overlap=True)

print("Loading low level intersections...", flush=True)
#CDS_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_CDS', bp_overlap=True)
#five_prime_UTR_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_five_prime_UTR', bp_overlap=True)
#three_prime_UTR_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_three_prime_UTR', bp_overlap=True)
#UTR_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, f'{ukb}/side_analyses/str_annotations/intersects_UTR', bp_overlap=True)
CDS_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, args.intersects_CDS_annotation, bp_overlap=True)
five_prime_UTR_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, args.intersects_five_prime_UTR_annotation, bp_overlap=True)
three_prime_UTR_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, args.intersects_three_prime_UTR_annotation, bp_overlap=True)
UTR_intersect_merge = annotation_utils.get_merged_annotations(finemapping_results, args.intersects_UTR_annotation, bp_overlap=True)

sub_exon_types = pd.concat([CDS_intersect_merge, five_prime_UTR_intersect_merge, three_prime_UTR_intersect_merge, UTR_intersect_merge])

nrows = finemapping_results.shape[0]

relation_to_gene = ['']*nrows

print("Handling gene and transcript annotations ...", flush=True)
for idx in range(nrows):
    chrom = finemapping_results['chrom'][idx]
    start_pos = finemapping_results['start_pos'][idx]
    end_pos = finemapping_results['end_pos'][idx]

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

    if len(gene_intersections) == 0:
        relation_to_gene[idx] = 'intergenic'

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

finemapping_results['relation_to_gene'] = relation_to_gene

# done with gene annotations, move back to polars
finemapping_results = pl.DataFrame(finemapping_results)

def apob_locus(phenotype):
    return (pl.col('start_pos') == 21266752) & (pl.col('phenotype') == phenotype)

finemapping_results = finemapping_results.with_columns([
    pl.when(~apob_locus('apolipoprotein_b')).then(pl.col('susie_CP')).otherwise(0.994568),
    pl.when(
        ~apob_locus('apolipoprotein_b') & ~apob_locus('ldl_cholesterol_direct')
    ).then(
        pl.col('susie_CP_best_guess_genotypes')
    ).when(apob_locus('apolipoprotein_b')).then(0.9765144).otherwise(0.9751951),
    pl.when(~apob_locus('apolipoprotein_b')).then(pl.col('finemap_CP_prior_4_signals')).otherwise(0.884233 + 0.0728702),
    pl.when(~apob_locus('apolipoprotein_b')).then(pl.col('susie_CP_prior_snps_over_strs')).otherwise(0.9779266),
    pl.when(~apob_locus('ldl_cholesterol_direct')).then(pl.col('finemap_CP_prior_snps_over_strs')).otherwise(1),
    pl.when(~apob_locus('apolipoprotein_b')).then(pl.col('finemap_CP_prior_effect_size_0.0025%')).otherwise(0.347805 + 0.652195),
])

finemapping_results = finemapping_results.select([
    'phenotype',
    'chrom',
    pl.col('start_pos').alias('start_pos (hg19)'),
    pl.col('end_pos').alias('end_pos (hg19)'),
    pl.col('pos_hg38').alias('start_pos (hg38)'),
    pl.col('end_pos_hg38').alias('end_pos (hg38)'),
    pl.when(~pl.col('pos_t2t').is_null()).then(pl.col('pos_t2t').cast(int).cast(str)).otherwise(pl.lit('Failed to lift to T2T reference')).alias('start_pos (T2T)'),
    pl.when(~pl.col('end_pos_t2t').is_null()).then(pl.col('end_pos_t2t').cast(int).cast(str)).otherwise(pl.lit('Failed to lift to T2T reference')).alias('end_pos (T2T)'),
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
        (pl.col('finemap_CP_stopping_thresh_1e-4') >= 0.8) &
        (pl.col('p_val') <= 1e-10)
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
    sum(
        [pl.col(f'{ethnicity}_p_val').cast(str) + pl.lit(', ') for ethnicity in other_ethnicities],
        pl.lit('')
    ).str.replace(', $', '').alias('other_ethnicity_association_p_values'),
    sum(
        [pl.when(pl.col(f'{ethnicity}_p_val') > .05).then('NA').when(pl.col('coeff') > 0).then('+').otherwise('-') + pl.lit(', ')
        for ethnicity in other_ethnicities],
        pl.lit('')
    ).str.replace(', $', '').alias('other_ethnicity_effect_directions'),
    *[pl.col(f'{ethnicity}_allele_dosages').apply(dosages_to_frequencies).alias(f'{ethnicity}_allele_frequencies') for ethnicity in other_ethnicities],
])

#finemapping_results.write_csv(f'{ukb}/post_finemapping/results/singly_finemapped_strs_for_paper.tab', sep='\t')
#finemapping_results.sort(['chrom', 'start_pos']).write_csv(f'{ukb}/post_finemapping/results/singly_finemapped_strs_sorted.tab', sep='\t')
finemapping_results.write_csv(f'{args.outdir}/singly_finemapped_strs_for_paper.tab', sep='\t')
finemapping_results.sort(['chrom', 'start_pos (hg19)']).write_csv(f'{args.outdir}/singly_finemapped_strs_sorted.tab', sep='\t')

confident_results = finemapping_results.filter(
    (pl.col('finemapping') == 'confidently').any().over(['chrom', 'start_pos (hg19)']) &
    (pl.col('association_p_value') <= 1e-10)
)
confident_results.write_csv(f'{args.outdir}/confidently_finemapped_strs_for_paper.tab', sep='\t')
confident_results.sort(['chrom', 'start_pos (hg19)']).write_csv(f'{args.outdir}/confidently_finemapped_strs_sorted.tab', sep='\t')
