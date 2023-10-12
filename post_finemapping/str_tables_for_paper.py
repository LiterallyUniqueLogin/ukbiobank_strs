#!/usr/bin/env python3

import argparse
import ast
import json
import os

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
parser.add_argument('--flank-start-to-start-and-end-pos')
parser.add_argument('--str-hg19-pos-bed')
parser.add_argument('--str-hg38-pos-bed')
parser.add_argument('--repeat-units-table')
parser.add_argument('--intersects-gene-annotation', nargs=22)
parser.add_argument('--intersects-exon-annotation', nargs=22)
parser.add_argument('--intersects-CDS-annotation', nargs=22)
parser.add_argument('--intersects-five-prime-UTR-annotation', nargs=22)
parser.add_argument('--intersects-three-prime-UTR-annotation', nargs=22)
parser.add_argument('--intersects-UTR-annotation', nargs=22)
parser.add_argument('--wgs-comparison-stats', nargs=6)
parser.add_argument('--wgs-allele-freqs', nargs=6)
parser.add_argument('--outdir')

args = parser.parse_args()

assert len(args.assoc_phenotypes) == \
        len(args.assocs) == \
        len(args.black_assocs) == \
        len(args.south_asian_assocs) == \
        len(args.chinese_assocs) == \
        len(args.irish_assocs) == \
        len(args.white_other_assocs)

flank_start_to_start_and_pos_table = pl.read_csv(
    args.flank_start_to_start_and_end_pos, sep='\t'
).groupby(
    ['chrom', 'pos', 'end_pos']
).agg([pl.col('snpstr_pos').first()]) # drop duplicate snpstr poses

hg19_pos_bed = pl.read_csv(
    args.str_hg19_pos_bed,
    sep='\t',
    has_header=False,
    new_columns=['chrom', 'pos', 'end_pos', 'ID']
).with_columns([
    pl.col('pos') + 1, # transform from bed to vcf coords
    pl.col('chrom').str.replace('chr', '').cast(int),
])
hg38_pos_bed = pl.read_csv(
    args.str_hg38_pos_bed,
    sep='\t',
    has_header=False,
    new_columns=['chrom', 'pos', 'end_pos', 'ID', 'drop']
).with_columns([
    pl.col('pos') + 1, # transform from bed to vcf coords
]).drop(['drop', 'chrom'])

hg19_pos_bed = hg19_pos_bed.join(
    hg38_pos_bed,
    how = 'left',
    on=['ID'],
    suffix='_hg38',
)

flank_start_to_start_and_pos_table = hg19_pos_bed.join(
    flank_start_to_start_and_pos_table,
    on=['chrom', 'pos', 'end_pos']
).drop(['ID'])


flank_start_to_start_and_pos_table = flank_start_to_start_and_pos_table.lazy()

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
for count, fname in enumerate(args.first_pass_finemapping_dfs):
    print(f'Loading finemapping df {count+1}/{len(args.first_pass_finemapping_dfs)}', flush=True)
    df = pl.scan_csv(
        fname,
        sep='\t',
        dtypes={
            **{f'{ethnicity}_p_val': float for ethnicity in other_ethnicities},
            **{f'{ethnicity}_coeff': float for ethnicity in other_ethnicities},
            **{f'{ethnicity}_se': float for ethnicity in other_ethnicities}
        }
    ).filter('is_STR')
    #phenotype = df['phenotype'][0] # all will be the same
    phenotype = fname.split('concordance_')[1].split('.tab')[0]
    df = df.filter(~(
    (
        (phenotype == 'total_bilirubin') &
        (pl.col('region') == '12_19976272_22524428')
    ) |
    (
        (phenotype == 'urate') &
        (pl.col('region') == '4_8165642_11717761')
    ) |
    (
        (phenotype == 'alkaline_phosphatase') &
        (pl.col('region') == '1_19430673_24309348')
    )))

    assoc_df = pl.scan_csv(
        assocs[phenotype],
        sep='\t',
    ).select([
        'chrom',
        'pos',
        pl.col('subset_total_per_allele_dosages').alias('white_british_allele_dosages')
    ])
    lazy_df = df.join(
        assoc_df,
        how='left',
        on=['chrom', 'pos']
    )
    for ethnicity in other_ethnicities:
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
    # hard code which phenotypes are using which coordinates
    if phenotype in {'ldl_cholesterol_direct', 'total_bilirubin'}:
        lazy_df = lazy_df.join(
            flank_start_to_start_and_pos_table,
            how='left',
            on=['chrom', 'pos']
        )
    else:
        lazy_df = lazy_df.rename({'pos': 'snpstr_pos'}).join(
            flank_start_to_start_and_pos_table,
            how='left',
            on=['chrom', 'snpstr_pos']
        )
    lazy_df = lazy_df.select([*[column for column in lazy_df.columns if column not in {'pos', 'snpstr_pos'}], 'pos', 'snpstr_pos'])
    finemapping_dfs.append(lazy_df)
finemapping_results = pl.concat(finemapping_dfs).collect()

finemapping_results = finemapping_results.filter(
    (pl.col('p_val') < 5e-8)
).filter(
    (
        ((pl.col('susie_alpha') >= 0.8) & (pl.col('susie_cs') >= 0)) | (pl.col('finemap_pip') >= 0.8)
    ).any().over(['chrom', 'snpstr_pos'])
)

repeat_units = pl.read_csv(
    args.repeat_units_table,
    sep='\t'
)

finemapping_results = finemapping_results.join(
    repeat_units,
    how='left',
    on=['chrom', 'snpstr_pos']
)

concordance_cols = pl.read_csv(
    args.followup_finemapping_dfs[0],
    sep='\t',
    n_rows=1
).columns

del lazy_df
followup_lazy_dfs = []
for count, followup_finemapping_df in enumerate(args.followup_finemapping_dfs):
    if os.stat(followup_finemapping_df).st_size <= 100:
        continue
    followup_lazy_df = pl.scan_csv(
        followup_finemapping_df,
        sep='\t',
        dtypes={col: (float if 'cs' not in col else int) for col in concordance_cols if 'finemap' in col or 'susie' in col or 'p_val' in col}
    ).filter('is_STR')

    # hard code which phenotypes are using which coordinates
    if phenotype in {'ldl_cholesterol_direct', 'total_bilirubin'}:
        followup_lazy_df = followup_lazy_df.join(
            flank_start_to_start_and_pos_table,
            how='left',
            on=['chrom', 'pos']
        )
    else:
        followup_lazy_df = followup_lazy_df.rename({'pos': 'snpstr_pos'}).join(
            flank_start_to_start_and_pos_table,
            how='left',
            on=['chrom', 'snpstr_pos']
        )
    followup_lazy_df = followup_lazy_df.select([*[column for column in followup_lazy_df.columns if column not in {'pos', 'snpstr_pos'}], 'pos', 'snpstr_pos'])
    followup_lazy_dfs.append(followup_lazy_df)
concordance_results = pl.concat(followup_lazy_dfs).collect()

finemapping_results = finemapping_results.join(
    concordance_results,
    how='left',
    on=['phenotype', 'chrom', 'snpstr_pos']
).with_columns([
    pl.when(pl.col('susie_alpha').is_null()).then(None).when(pl.col('susie_cs') >= 0).then(pl.col('susie_alpha')).otherwise(0).alias('susie_CP'),
    pl.when(pl.col('susie_alpha_best_guess').is_null()).then(None).when(pl.col('susie_cs_best_guess') >= 0).then(pl.col('susie_alpha_best_guess')).otherwise(0).alias('susie_CP_best_guess_genotypes'),
    pl.when(pl.col('susie_alpha_ratio').is_null()).then(None).when(pl.col('susie_cs_ratio') >= 0).then(pl.col('susie_alpha_ratio')).otherwise(0).alias('susie_CP_prior_snps_over_strs'),
    pl.col('finemap_pip').alias('finemap_CP'),
    pl.col('finemap_pip_repeat').alias('finemap_CP_repeat'),
    pl.col('finemap_pip_p_thresh').alias('finemap_CP_pval_thresh_5e-4'),
    pl.col('finemap_pip_mac').alias('finemap_CP_mac_thresh_100'),
    pl.col('finemap_pip_prior_std_derived').alias('finemap_CP_prior_effect_size_0.05%'),
    pl.col('finemap_pip_total_prob').alias('finemap_CP_prior_4_signals'),
    pl.col('finemap_pip_conv_tol').alias('finemap_CP_stopping_thresh_1e-4'),
    pl.col('finemap_pip_ratio').alias('finemap_CP_prior_snps_over_strs'),
    pl.col('finemap_pip_prior_std_low').alias('finemap_CP_prior_effect_size_0.0025%'),
    pl.col('pos').alias('start_pos')
])

# WGS
def clean_allele_freqs(freq_dict_str):
    allele_freqs = ast.literal_eval('{' + freq_dict_str + '}')
    return json.dumps({ (int(k) if k == int(k) else round(k, 2)) : f'{v*100:.1f}%' for k,v in allele_freqs.items() }).replace('"', '')

def make_len_sum_column(cols):
    len_sum_frequencies = ast.literal_eval(cols['len_sum_frequencies'])
    len_sum_accuracies = ast.literal_eval(cols['len_sum_accuracies'])
    assert set(len_sum_frequencies) == set(len_sum_accuracies)
    return json.dumps({ int(k): f'({len_sum_frequencies[k]*100:.2f}%, {len_sum_accuracies[k]*100:.2f}%)' for k in len_sum_frequencies}).replace('"', '')

for ethnicity, allele_freqs_fname, comparison_stats_fname in zip(
    ['white_british'] + other_ethnicities,
    args.wgs_allele_freqs,
    args.wgs_comparison_stats
):
    freqs = pl.read_csv(
        allele_freqs_fname,
        sep='\t'
    ).select([
        pl.col('afreq-1').apply(clean_allele_freqs).alias(f'{ethnicity}_WGS_allele_frequencies'),
        pl.col('chrom').str.replace('chr', '').cast(int),
        pl.col('start').alias('pos_hg38'),
    ])
    finemapping_results = finemapping_results.join(freqs, on=['chrom', 'pos_hg38'])

    comparison_stats = pl.read_csv(
        comparison_stats_fname,
        sep='\t'
    ).select([
        'chrom',
        pl.col('start').alias('start_pos'),
        pl.col('numcalls').alias(f'{ethnicity}_WGS_numcalls'),
        pl.col('fraction_concordant_len_sum').round(4).alias(f'{ethnicity}_fraction_concordant_len_sum'),
        pl.col('mean_absolute_difference').round(4).alias(f'{ethnicity}_mean_absolute_difference'),
        pl.col('r').round(4).alias(f'{ethnicity}_r'),
        pl.col('dosage_r').round(4).alias(f'{ethnicity}_dosage_r'),
        pl.struct(['len_sum_frequencies', 'len_sum_accuracies']).apply(make_len_sum_column).alias(f'{ethnicity}_length_sum_accuracies')
    ])
    finemapping_results = finemapping_results.join(comparison_stats, on=['chrom', 'start_pos'])

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

finemapping_results = finemapping_results.select([
    'phenotype',
    'chrom',
    pl.col('start_pos').alias('start_pos (hg19)'),
    pl.col('end_pos').alias('end_pos (hg19)'),
    pl.col('pos_hg38').alias('start_pos (hg38)'),
    pl.col('end_pos_hg38').alias('end_pos (hg38)'),
    pl.col('region').alias('finemapping_region'),
    pl.when(pl.col('unit') != 'None').then(pl.col('unit')).otherwise('None:' + pl.col('period').cast(str)).alias('repeat_unit'),
    pl.col('white_british_allele_dosages').apply(dosages_to_frequencies).alias('white_british_imputed_allele_frequencies'),
    pl.col('p_val').alias('association_p_value'),
    pl.when(pl.col('coeff') > 0).then('+').otherwise('-').alias('direction_of_association'),
    'relation_to_gene',
    pl.when(
        (pl.col('susie_CP') >= 0.8) &
        (pl.col('finemap_CP') >= 0.8) &
        (pl.col('susie_CP_best_guess_genotypes') >= 0.8) &
        (pl.col('finemap_CP_repeat') >= 0.8) &
        (pl.col('finemap_CP_pval_thresh_5e-4') >= 0.8) &
        (pl.col('finemap_CP_mac_thresh_100') >= 0.8) &
        (pl.col('finemap_CP_prior_effect_size_0.05%') >= 0.8) &
        (pl.col('finemap_CP_prior_4_signals') >= 0.8) &
        (pl.col('finemap_CP_stopping_thresh_1e-4') >= 0.8) &
        (pl.col('p_val') < 1e-10)
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
    'finemap_CP_repeat',
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
    *[pl.col(f'{ethnicity}_allele_dosages').apply(dosages_to_frequencies).alias(f'{ethnicity}_imputed_allele_frequencies') for ethnicity in other_ethnicities],
    sum(
        [pl.col(f'{ethnicity}_WGS_numcalls') + pl.lit(', ') for ethnicity in ['white_british'] + other_ethnicities],
        pl.lit('')
    ).str.replace(', $', '').alias('num_WGS_calls_per_population'),
    *[f'{ethnicity}_WGS_allele_frequencies' for ethnicity in ['white_british'] + other_ethnicities],
    sum(
        [pl.col(f'{ethnicity}_fraction_concordant_len_sum') + pl.lit(', ') for ethnicity in ['white_british'] + other_ethnicities],
        pl.lit('')
    ).str.replace(', $', '').alias('fraction_WGS_concordant_length_sums'),
    sum(
        [pl.col(f'{ethnicity}_mean_absolute_difference') + pl.lit(', ') for ethnicity in ['white_british'] + other_ethnicities],
        pl.lit('')
    ).str.replace(', $', '').alias('mean_absolute_difference_with_WGS'),
    sum(
        [pl.col(f'{ethnicity}_r') + pl.lit(', ') for ethnicity in ['white_british'] + other_ethnicities],
        pl.lit('')
    ).str.replace(', $', '').alias('WGS_weighted_r'),
    sum(
        [pl.col(f'{ethnicity}_dosage_r') + pl.lit(', ') for ethnicity in ['white_british'] + other_ethnicities],
        pl.lit('')
    ).str.replace(', $', '').alias('WGS_dosage_r'),
    *[pl.col(f'{ethnicity}_length_sum_accuracies').alias(f'{ethnicity}_WGS_per_length_sum_concordances') for ethnicity in ['white_british'] + other_ethnicities]
])

finemapping_results.write_csv(f'{args.outdir}/singly_finemapped_strs_for_paper.tab', sep='\t')
finemapping_results.sort(['chrom', 'start_pos (hg19)']).write_csv(f'{args.outdir}/singly_finemapped_strs_sorted.tab', sep='\t')

confident_results = finemapping_results.filter(
    (pl.col('association_p_value') < 1e-10)
).filter(
    (pl.col('finemapping') == 'confidently').any().over(['chrom', 'start_pos (hg19)'])
)
confident_results.write_csv(f'{args.outdir}/confidently_finemapped_strs_for_paper.tab', sep='\t')
confident_results.sort(['chrom', 'start_pos (hg19)']).write_csv(f'{args.outdir}/confidently_finemapped_strs_sorted.tab', sep='\t')
