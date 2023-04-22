#!/usr/bin/env python3

import argparse

import numpy as np
import polars as pl

import annotation_utils

other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']

parser = argparse.ArgumentParser()
parser.add_argument('--outdir')
parser.add_argument('--str-pos-table')
parser.add_argument('--str-loci')
parser.add_argument('--repeat-units')
parser.add_argument('--eSTR-table')
parser.add_argument('--gencode')
parser.add_argument('--phenotypes', nargs='+')
parser.add_argument('--assocs', nargs='+')
parser.add_argument('--confidently-finemapped-STRs-df')
parser.add_argument('--intersects-protein-coding-CDS-support-2', nargs=22)
parser.add_argument('--intersects-protein-coding-UTR-support-2', nargs=22)
parser.add_argument('--intersects-protein-coding-five-prime-UTR-support-2', nargs=22)
parser.add_argument('--intersects-protein-coding-three-prime-UTR-support-2', nargs=22)
parser.add_argument('--intersects-transcript-support-2', nargs=22)
parser.add_argument('--closest-downstream-protein-coding-exon-support-2', nargs=22)
parser.add_argument('--closest-upstream-protein-coding-exon-support-2', nargs=22)
parser.add_argument('--closest-downstream-protein-coding-gene', nargs=22)
parser.add_argument('--closest-upstream-protein-coding-gene', nargs=22)
# TODO causal STR candidates - add filter?
args = parser.parse_args()

# pos (start), snpstr_pos (hipstr)
all_STRs = pl.read_csv(
    #f'{ukb}/snpstr/flank_trimmed_vcf/vars.tab',
    args.str_pos_table,
    sep='\t'
)
# pos (hisptr)
snpstr_strs = pl.scan_csv(
    #f'{ukb}/snpstr/str_loci.txt',
    args.str_loci,
    sep='\t',
    has_header=False,
    with_column_names = lambda _ : ['chrom', 'pos'],
)

all_STRs = all_STRs.lazy().join(
    snpstr_strs,
    left_on=['chrom', 'snpstr_pos'],
    right_on=['chrom', 'pos'],
    how='inner',
    suffix='_other'
).select(
    ['chrom', 'pos', 'end_pos', 'snpstr_pos']
).with_column(
    pl.col('snpstr_pos').alias('SNPSTR_start_pos')
).drop('snpstr_pos').unique(
    subset=['chrom', 'pos']
).collect()
assert ~np.any(np.isnan(all_STRs['chrom'].to_numpy()))
assert ~np.any(np.isnan(all_STRs['pos'].to_numpy()))
assert ~np.any(np.isnan(all_STRs['end_pos'].to_numpy()))
assert ~np.any(np.isnan(all_STRs['SNPSTR_start_pos'].to_numpy()))
# pos (start), end_pos, SNPSTR_start_pos (hipstr)

repeat_units = pl.read_csv(
    #f'{ukb}/snpstr/repeat_units.tab',
    args.repeat_units,
    sep='\t',
    columns=['chrom', 'snpstr_pos', 'period', 'unit']
).with_column(
    pl.col('unit').alias('canonical_unit')
).drop('unit')

all_STRs = all_STRs.join(
    repeat_units,
    left_on=['chrom', 'SNPSTR_start_pos'],
    right_on=['chrom', 'snpstr_pos'],
    how='left',
    suffix='_other'
)
assert ~np.any(np.isnan(all_STRs['period'].to_numpy()))

print('Loading eSTRs ... ')
eSTRs = pl.read_csv(
    #f'{ukb}/misc_data/eSTR/eSTRs.csv',
    args.eSTR_table,
    sep=','
).rename({
    'score': 'eSTR_CAVIAR_score'
}).with_column(
    pl.col('chrom').str.slice(3).cast(int)
).groupby(
    ['chrom', 'str.start']
).agg(pl.col('eSTR_CAVIAR_score').max())

all_STRs = all_STRs.join(
    eSTRs,
    how='left',
    left_on=['chrom', 'pos'],
    right_on=['chrom', 'str.start'],
).with_columns([
    (~pl.col('eSTR_CAVIAR_score').is_null()).alias('eSTR'),
    pl.when(
        pl.col('eSTR_CAVIAR_score').is_null()
    ).then(
        False
    ).otherwise(
        pl.col('eSTR_CAVIAR_score') >= .3
    ).alias('FM_eSTR')
])

print('Getting promoters ... ', flush=True, end='')
genes = pl.read_csv(
    #f'{ukb}/misc_data/gencode/gencode.v38lift37.annotation.without_chr.sorted.gene.gff3',
    args.gencode,
    sep='\t',
    has_header = False,
    columns=[0,3,4,6,8],
    dtypes={'column_1': str}
).select([
    pl.col('column_1').alias('chrom'),
    pl.col('column_4').alias('start_pos'),
    pl.col('column_5').alias('end_pos'),
    pl.col('column_7').alias('strand'),
    pl.col('column_9').alias('kvps')
]).with_column(
    pl.when(
        pl.col('strand') == '+'
    ).then(
        pl.col('start_pos')
    ).when(
        pl.col('strand') == '-'
    ).then(
        pl.col('end_pos')
    ).otherwise(
        pl.lit(None)
    ).alias('TSS')
).filter(
    pl.col('chrom').is_in([str(x) for x in range(1,23)])
)
assert not np.any(np.isnan(genes['TSS'].to_numpy()))

all_STRs = all_STRs.to_pandas()
all_STRs['upstream_promoter'] = False
for row in range(genes.shape[0]):
    if not 'gene_type=protein_coding' in genes[row, 'kvps']:
        continue
    if genes[row, 'chrom'] not in set(str(x) for x in range(1, 23)):
        continue
    tss = genes[row, 'TSS']
    all_STRs.loc[
        (all_STRs['chrom'] == int(genes[row, 'chrom'])) & (
            ((tss - all_STRs['end_pos'] <= 3000) & (tss - all_STRs['end_pos'] >= 0) & (genes[row, 'strand'] == '+')) |
            ((all_STRs['pos'] - tss <= 3000) & (all_STRs['pos'] - tss >= 0) & (genes[row, 'strand'] == '-'))
        ),
        'upstream_promoter'
    ] = True
print('done', flush=True)

#annotation_dir = f'{ukb}/side_analyses/str_annotations'

for fnames, annotation_type in [
    (args.intersects_protein_coding_CDS_support_2, 'coding'),
    (args.intersects_protein_coding_UTR_support_2, 'UTR'),
    (args.intersects_protein_coding_five_prime_UTR_support_2, 'UTR5'),
    (args.intersects_protein_coding_three_prime_UTR_support_2, 'UTR3'),
    (args.intersects_transcript_support_2, 'transcribed')
]:
    print(f'Loading annotation {annotation_type} ... ', flush=True, end='')
    intersects_df = annotation_utils.get_merged_annotations(
        all_STRs, fnames, how='left', bp_overlap=True
    ).drop_duplicates(subset=['chrom', 'pos']).reset_index(drop=True)
    # pos (start), SNPSTR_start_pos (hipstr)
    # could intersect multiple things, so need to group by to see if intersecting any
    # https://stackoverflow.com/questions/43321455/pandas-count-null-values-in-a-groupby-function
    # start with any field from the annotation df
    intersects_df = intersects_df['annotation_pos'] \
        .isnull() \
        .groupby([intersects_df['chrom'], intersects_df['pos']]) \
        .sum() \
        .astype(bool) \
        .reset_index()
    assert np.all(intersects_df[['chrom', 'pos']] == all_STRs[['chrom', 'pos']])
    all_STRs[annotation_type] = ~intersects_df['annotation_pos']
    print('done', flush=True)

for fnames, annotation_type in [
    (args.closest_downstream_protein_coding_exon_support_2, 'downstream_exon'),
    (args.closest_upstream_protein_coding_exon_support_2, 'upstream_exon'),
    (args.closest_downstream_protein_coding_gene, 'downstream_gene'),
    (args.closest_upstream_protein_coding_gene, 'upstream_gene')
]:
    print(f'Loading annotation {annotation_type} ... ', flush=True, end='')
    closest_df = annotation_utils.get_merged_annotations(
        all_STRs, fnames, how='left', distance=True
    ).drop_duplicates(subset=['chrom', 'pos']).reset_index(drop=True)
    assert np.all(closest_df[['chrom', 'pos']] == all_STRs[['chrom', 'pos']])
    all_STRs[f'{annotation_type}_dist'] = np.abs(closest_df['annotation_distance'])
    print('done', flush=True)

'''
#intersection should be the same regardless of direction
assert np.all((all_STRs['downstream_exon_dist'] == 0) == (all_STRs['upstream_exon_dist'] == 0))
assert np.all((all_STRs['downstream_gene_dist'] == 0) == (all_STRs['upstream_gene_dist'] == 0))
'''

all_STRs['exonic'] = all_STRs['downstream_exon_dist'] == 0
all_STRs['genic'] = all_STRs['downstream_gene_dist'] == 0
for key in 'exonic', 'genic':#, 'UTR', 'UTR3', 'UTR5':#, 'promoter':
    print(f'n {key} STRS: {np.sum(all_STRs[key])}')
    print(f'n non-genic {key} STRS: {np.sum(~all_STRs.loc[all_STRs[key], "genic"])}')

# not all exonic STRs should be coding - UTRs are also considered exons
# all UTRs and exons should be genic
assert np.all(all_STRs.loc[all_STRs['exonic'], 'genic'])
assert np.all(all_STRs.loc[all_STRs['UTR'], 'genic'])
assert np.all(all_STRs.loc[all_STRs['UTR3'], 'genic'])
assert np.all(all_STRs.loc[all_STRs['UTR5'], 'genic'])
# also exlcude 23 introns that aren't marked as exons or UTRs but are marked as coding
all_STRs['intronic'] = (
    all_STRs['genic'] & ~all_STRs['exonic'] & ~all_STRs['UTR'] & ~all_STRs['UTR5'] & ~all_STRs['UTR3'] & ~all_STRs['coding']
)
all_STRs['transcribed_non_protein'] = ~all_STRs['genic'] & all_STRs['transcribed']
all_STRs['intergenic'] = ~all_STRs['genic']
for key in ('intronic', 'transcribed_non_protein'):
    print(f'n {key} STRS: {np.sum(all_STRs[key])}')

causal_STR_candidates = pl.read_csv(
    args.confidently_finemapped_STRs_df,
    #f'{ukb}/post_finemapping/intermediate_results/concordant_causal_STR_candidates.tab',
    sep='\t'
).select([
    'phenotype',
    'chrom',
    'pos',
    pl.lit(True).alias('is_causal_STR_candidate')
])

all_STRs = pl.DataFrame(all_STRs).join(
    causal_STR_candidates.groupby(['chrom', 'pos']).agg(pl.col('is_causal_STR_candidate').any()),
    how='left',
    left_on=['chrom', 'SNPSTR_start_pos'],
    right_on=['chrom', 'pos']
)

assert all_STRs.groupby(['chrom', 'pos']).agg(pl.count()).select(pl.col('count').max().alias('out'))['out'].to_numpy()[0] == 1

'''
finemapping_results = pl.read_csv(
    'post_finemapping/intermediate_results/finemapping_putatively_causal_concordance.tab',
    sep='\t'
).filter(
    ~pl.col('finemap_pip').is_null() &
    ~pl.col('susie_alpha').is_null() &
    pl.col('is_STR') &
    (pl.col('p_val') <= 1e-10)
).with_columns([
    pl.when(pl.col('susie_cs') > 0).then(pl.col('susie_alpha')).otherwise(0).alias('susie_alpha'),
    pl.when(pl.col('susie_cs_ratio') > 0).then(pl.col('susie_alpha_ratio')).otherwise(0).alias('susie_alpha_ratio'),
    pl.when(pl.col('susie_cs_hardcall') > 0).then(pl.col('susie_alpha_hardcall')).otherwise(0).alias('susie_alpha_hardcall'),
])

susie_cols = finemapping_results.select([
        pl.col('^susie_alpha.*$'),
]).columns
finemap_cols = finemapping_results.select([
        pl.col('^finemap_pip.*$')
]).columns

either = finemapping_results.filter(
    (pl.col('susie_alpha') >= .8) |
    (pl.col('finemap_pip') >= .8)
).select(['chrom', 'pos']).distinct().with_column(pl.lit(True).alias('indicator'))

both = finemapping_results.filter(
    (pl.col('susie_alpha') >= .8) &
    (pl.col('finemap_pip') >= .8)
).select(['chrom', 'pos']).distinct().with_column(pl.lit(True).alias('indicator'))

resilient_but_ratio_low = finemapping_results.filter(
    (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols if 'ratio' not in col and 'prior_std_low' not in col]) == 8)
).select(['chrom', 'pos']).distinct().with_column(pl.lit(True).alias('indicator'))
print(
    'n STR-phenotype assocs',
    finemapping_results.filter(
        (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols if 'ratio' not in col and 'prior_std_low' not in col]) == 8)
    ).shape
)

resilient = finemapping_results.filter(
    (pl.sum([(pl.col(col) >= .8).cast(int) for col in susie_cols + finemap_cols]) == 11)
).select(['chrom', 'pos']).distinct().with_column(pl.lit(True).alias('indicator'))
#print(either.shape[0], both.shape[0], resilient_but_ratio_low.shape[0], resilient.shape[0])

# subset to STRs that were fine-mapped by both fine-mappers in a pheno with p-val <= 1e-10, then
# mark which ones passed different fine-mapping thresholds
compare_STRs = pl.DataFrame(all_STRs).join(
    finemapping_results.select(['chrom', 'pos']).distinct(),
    how='inner',
    left_on = ['chrom', 'SNPSTR_start_pos'],
    right_on = ['chrom', 'pos'],
).join(
    either,
    how='left',
    left_on = ['chrom', 'SNPSTR_start_pos'],
    right_on = ['chrom', 'pos'],
).with_column(
    (~pl.col('indicator').is_null()).alias('singly_finemapped')
).drop('indicator').join(
    both,
    how='left',
    left_on = ['chrom', 'SNPSTR_start_pos'],
    right_on = ['chrom', 'pos'],
).with_column(
    (~pl.col('indicator').is_null()).alias('doubly_finemapped')
).drop('indicator').join(
    resilient_but_ratio_low,
    how='left',
    left_on = ['chrom', 'SNPSTR_start_pos'],
    right_on = ['chrom', 'pos'],
).with_column(
    (~pl.col('indicator').is_null()).alias('concordantly_finemapped')
).drop('indicator').join(
    resilient,
    how='left',
    left_on = ['chrom', 'SNPSTR_start_pos'],
    right_on = ['chrom', 'pos'],
).with_column(
    (~pl.col('indicator').is_null()).alias('strictly_concordantly_finemapped')
).drop('indicator')
'''

print('loading p_vals ... ', flush=True)
gwsig_loci = set()
for phenotype, fname in zip(args.phenotypes, args.assocs):
    print(phenotype, flush=True)
    #fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    with open(fname) as tsv:
        header = tsv.readline().strip()
    gwsig_df = pl.scan_csv(
        fname,
        sep = '\t',
        skip_rows=1,
        has_header=False,
        with_column_names = lambda _: header.replace('0.05_significance_CI', 'foo', 1).replace('5e-8_significance_CI', 'bar', 1).split('\t') # these duplicate column names won't be used anyway
    ).filter(pl.col(f'p_{phenotype}') < 5e-8).select(['chrom', 'pos']).collect()
    gwsig_loci = gwsig_loci.union(zip(gwsig_df['chrom'].to_numpy(), gwsig_df['pos'].to_numpy()))
gwsig_loci = [f'{chrom}_{pos}' for (chrom, pos) in gwsig_loci]
print('matching gwsig loci .... ', flush=True)
all_STRs = all_STRs.with_column(
    (pl.col('chrom').cast(str) + '_' +  pl.col('SNPSTR_start_pos').cast(str)).is_in(gwsig_loci).alias('gwsig')
)

all_STRs.write_csv(
    f'{args.outdir}/enrichments_df.tab',
    #f'{ukb}/post_finemapping/intermediate_results/enrichments_df.tab',
    sep='\t'
)
