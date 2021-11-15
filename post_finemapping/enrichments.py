#!/usr/bin/env python3

import argparse
import ast
import os
from typing import List

import bokeh.io
import bokeh.models
import bokeh.models.tickers
import bokeh.plotting
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import scipy.stats.contingency

import annotation_utils
import str_utils

ukb = os.environ['UKB']

'''
ALL STRs
And with SIG <= 5e-8
And with FINEMAP pcausal >= .9
'''
#def main():

def binary_test(out, compare_STRs, in_category, category_string):
    FINEMAPed = compare_STRs['FINEMAPed']
    # marginals
    n_strs = compare_STRs.shape[0]
    n_cat = np.sum(in_category)
    n_subset = np.sum(FINEMAPed)
    # table elements
    n_neither = np.sum(~FINEMAPed & ~in_category)
    n_subset_not_cat = np.sum(FINEMAPed & ~in_category)
    n_cat_not_subset = np.sum(~FINEMAPed & in_category)
    n_subset_cat = np.sum(FINEMAPed & in_category)
    contingency_table = [[n_neither, n_subset_not_cat], [n_cat_not_subset, n_subset_cat]]
    if np.any(np.array(contingency_table) < 20):
        p_val = scipy.stats.fisher_exact(contingency_table)[1]
    else:
        p_val = scipy.stats.chi2_contingency(contingency_table)[1]
    out.write(f'{category_string}\t\t{n_strs}\t{n_cat}\t{n_cat/n_strs*100:.4f}%\t{n_subset}\t{n_subset_cat}\t{n_subset_cat/n_subset*100:.4f}%\t{p_val}\n')

def graph_relative_rate_by_pcausal(compare_STRs, in_category, category_string):
    xs = np.arange(0.00, 1.01, 0.01)
    finemaped_counts = []
    finemaped_subset_counts = []
    comparison_count = compare_STRs.shape[0]
    comparison_subset_count = np.sum(in_category)
    rate_ratios = []
    upper_cis = []
    lower_cis = []

    for x in xs:
        finemaped_count = np.sum(compare_STRs['FINEMAP_pcausal'] >= x)
        finemaped_counts.append(finemaped_count)
        finemaped_subset_count = np.sum(in_category & (compare_STRs['FINEMAP_pcausal'] >= x))
        finemaped_subset_counts.append(finemaped_subset_count)
        result = scipy.stats.contingency.relative_risk(
            finemaped_subset_count,
            finemaped_count,
            comparison_subset_count,
            comparison_count
        )
        rate_ratios.append(result.relative_risk)
        # +/-1 SD
        lower, upper = result.confidence_interval(.68)
        lower_cis.append(lower)
        upper_cis.append(upper)
    rate_ratios = np.array(rate_ratios)
    finemaped_counts = np.array(finemaped_counts)
    finemaped_subset_counts = np.array(finemaped_subset_counts)

    fig = bokeh.plotting.figure(
        title=f'Relative rate of {category_string} STRs by FINEMAP posterior causality',
        y_axis_label = 'Relative rate',
        x_axis_label = 'FINEMAP posterior causality at least this',
        width=1200,
        height=800
    )
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.grid.grid_line_color = None
    fig.toolbar_location = None
    fig.title.text_font_size = '30px'
    fig.axis.axis_label_text_font_size = '26px'
    fig.axis.major_label_text_font_size = '20px'
    fig.varea(xs, lower_cis, upper_cis, legend_label='+/- 1 SD')
    fig.line(xs, rate_ratios, color='black', legend_label='relative rate')
    fig.legend.label_text_font_size = '22px'
    fig.add_layout(bokeh.models.Title(
        text=f'Relative rate is compared to the set of all non-MHC gwas-sig STRs: {comparison_subset_count}/{comparison_count} = {comparison_subset_count/comparison_count*100:.4f}%'
    ), 'below')
    fig.add_layout(bokeh.models.Title(
        text=f'Relative rate is defined as (% {category_string} STRs among FINEMAP STRs with posterior causality >= x) / (% {category_string} STRs among all non-MHC gwas-sig STRs)'
    ), 'below')
    fig.add_layout(bokeh.models.Title(
        text='gwas-sig STRs are those with association p <= 5e-8, not in the MHC. FINEMAP STRs must also be gwas-sig'
    ), 'below')
    label_idxs = np.array([0,20,40,60,80,100])
    fig.circle(x=xs[label_idxs], y=rate_ratios[label_idxs], size=10, color='black')
    fig.add_layout(bokeh.models.LabelSet(
        x='x', y='y', text='text', x_offset=.03,  y_offset=.03, text_color='black', source=bokeh.models.ColumnDataSource(data=dict(
            x=xs[label_idxs], y=rate_ratios[label_idxs], text=[
                f'{subset}/{total}' for (subset, total) in zip(finemaped_subset_counts[label_idxs], finemaped_counts[label_idxs])
            ]
        ))
    ))
    fig.line([0,1], [1,1], line_dash='dashed', line_width=2, color='black')
    bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/FINEMAP_{category_string}_relative_rate.png')
    bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/FINEMAP_{category_string}_relative_rate.svg')

def _BetterCDF(data_list: List[float],
           ax: matplotlib.axes.Axes):
    # assumes that axes are already set to (min, max)
    data = np.sort(data_list)
    x_axis_min, x_axis_max = ax.get_xlim()
    n_points = len(data)
    data = np.hstack((
        [x_axis_min],
        data,
        [x_axis_max]
    ))
    ys = 1 - np.hstack((
        [1],
        np.arange(n_points - 1, -1, -1) / n_points,
        [0]
    ))
    return ax.step(data, ys, where='post')

parser = argparse.ArgumentParser()
parser.add_argument('phenotypes', nargs='+')
args = parser.parse_args()
phenotypes = args.phenotypes

# pos (start), snpstr_pos (hipstr)
all_STRs = pd.read_csv(
    f'{ukb}/snpstr/flank_trimmed_vcf/vars.tab',
    header=0,
    delimiter='\t'
)
# pos (hisptr)
snpstr_strs = pd.read_csv(
    f'{ukb}/snpstr/str_loci.txt',
    header=None,
    names=['chrom', 'pos'],
    delimiter='\t'
)

all_STRs = all_STRs.merge(
    snpstr_strs,
    left_on=['chrom', 'snpstr_pos'],
    right_on=['chrom', 'pos'],
    how='inner',
    suffixes=['', '_other']
)
all_STRs = all_STRs[['chrom', 'pos', 'end_pos', 'snpstr_pos']]
all_STRs.rename(columns={'snpstr_pos': 'SNPSTR_start_pos'}, inplace=True)
all_STRs = all_STRs.drop_duplicates(subset=['chrom', 'pos']).reset_index(drop=True)
# pos (start), end_pos, SNPSTR_start_pos (hipstr)

repeat_units = pd.read_csv(
    f'{ukb}/snpstr/repeat_units.tab',
    delimiter='\t',
    header=0,
    usecols=['chrom', 'snpstr_pos', 'period', 'unit']
)
all_STRs = all_STRs.merge(
    repeat_units,
    left_on=['chrom', 'SNPSTR_start_pos'],
    right_on=['chrom', 'snpstr_pos'],
    how='left',
    suffixes=['', 'other']
)
assert ~np.any(all_STRs['period'].isnull())
all_STRs['canonical_unit'] = [str_utils.canonicalize(seq) if seq != 'None' else 'None' for seq in all_STRs['unit']]

dna_structures = pd.read_csv(
    f'{ukb}/misc_data/stalling/canon_structure.tab',
    delimiter='\t',
    index_col=None
)

all_STRs = all_STRs.merge(
    dna_structures,
    how='left',
    left_on=['canonical_unit'],
    right_on=['repeat_unit'],
    suffixes=['', '_other']
)

loci_summary_dfs = []
for chrom in range(1, 23):
    distribution_stats = pd.read_csv(
        f'{ukb}/export_scripts/intermediate_results/chr{chrom}_loci_summary.tab',
        delimiter='\t',
        header=0,
        index_col=None
    )
    loci_summary_dfs.append(distribution_stats)
loci_summaries = pd.concat(loci_summary_dfs)
n_before = all_STRs.shape[0]
all_STRs = all_STRs.merge(
    loci_summaries,
    how='left',
    left_on=['chrom', 'snpstr_pos'],
    right_on=['chr', 'pos'],
    suffixes=['', '_other']
)
assert n_before == all_STRs.shape[0]
print('Calculating mean lens ... ', flush=True, end='')
all_STRs['mean_len'] = [
    sum(key*val for (key, val) in ast.literal_eval(allele_dist).items()) for allele_dist in all_STRs['allele_dist']
]
print('done', flush=True)

eSTRs = pd.read_csv(
    f'{ukb}/misc_data/eSTR/eSTRs.csv',
    delimiter=',',
    index_col=None,
    header=0
)
eSTRs.rename(columns={'score': 'eSTR_CAVIAR_score'}, inplace=True)
eSTRs = eSTRs.groupby(['chrom', 'str.start']).max().reset_index()
# trim off thet leading chs
eSTRs['chrom'] = [
    int(chrom[3:]) for chrom in eSTRs['chrom']
]
all_STRs = all_STRs.merge(
    eSTRs,
    how='left',
    left_on=['chrom', 'pos'],
    right_on=['chrom', 'str.start'],
    suffixes=['', '_other']
)
all_STRs['eSTR'] = ~all_STRs['eSTR_CAVIAR_score'].isnull()
all_STRs['FM_eSTR'] = all_STRs['eSTR_CAVIAR_score'] > .3

print('Getting promoters ... ', flush=True, end='')
all_STRs['promoter'] = False
genes = pd.read_csv(
    f'{ukb}/misc_data/gencode/gencode.v38lift37.annotation.without_chr.sorted.gene.gff3',
    delimiter='\t',
    names=['chrom', 'start_pos', 'end_pos', 'strand', 'kvps'],
    usecols=[0,3,4,6,8],
    index_col=None
)
genes['TSS'] = np.nan
genes.loc[genes['strand'] == '+', 'TSS'] = genes.loc[genes['strand'] == '+', 'start_pos']
genes.loc[genes['strand'] == '-', 'TSS'] = genes.loc[genes['strand'] == '-', 'end_pos']
assert not np.any(genes['TSS'].isnull())
for row in range(genes.shape[0]):
    if not 'gene_type=protein_coding' in genes.loc[row, 'kvps']:
        continue
    if genes.loc[row, 'chrom'] not in set(str(x) for x in range(1, 23)):
        continue
    tss = genes.loc[row, 'TSS']
    all_STRs.loc[
        (all_STRs['chrom'] == int(genes.loc[row, 'chrom'])) & (
            (np.abs(all_STRs['pos']-tss) <= 3000) | (np.abs(all_STRs['end_pos']-tss) <= 3000)
        ),
        'promoter'
    ] = True
print('done', flush=True)

annotation_dir = f'{ukb}/side_analyses/str_annotations'

for fname, annotation_type in [
    ('intersects_protein_coding_CDS_support_2', 'coding'),
    ('intersects_protein_coding_UTR_support_2', 'UTR'),
    ('intersects_protein_coding_five_prime_UTR_support_2', 'UTR5'),
    ('intersects_protein_coding_three_prime_UTR_support_2', 'UTR3'),
    ('intersects_transcript_support_2', 'transcribed')
]:
    print(f'Loading annotation {annotation_type} ... ', flush=True, end='')
    intersects_df = annotation_utils.get_merged_annotations(
        all_STRs, f'{annotation_dir}/{fname}', how='left', bp_overlap=True
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

for fname, annotation_type in [
    ('closest_downstream_protein_coding_exon_support_2', 'downstream_exon'),
    ('closest_upstream_protein_coding_exon_support_2', 'upstream_exon'),
    ('closest_downstream_protein_coding_gene', 'downstream_gene'),
    ('closest_upstream_protein_coding_gene', 'upstream_gene')
]:
    print(f'Loading annotation {annotation_type} ... ', flush=True, end='')
    closest_df = annotation_utils.get_merged_annotations(
        all_STRs, f'{annotation_dir}/{fname}', how='left', distance=True
    ).drop_duplicates(subset=['chrom', 'pos']).reset_index(drop=True)
    assert np.all(closest_df[['chrom', 'pos']] == all_STRs[['chrom', 'pos']])
    all_STRs[f'{annotation_type}_dist'] = np.abs(closest_df['annotation_distance'])
    print('done', flush=True)

#intersection should be the same regardless of direction
assert np.all((all_STRs['downstream_exon_dist'] == 0) == (all_STRs['upstream_exon_dist'] == 0))
assert np.all((all_STRs['downstream_gene_dist'] == 0) == (all_STRs['upstream_gene_dist'] == 0))

all_STRs['exonic'] = all_STRs['downstream_exon_dist'] == 0
all_STRs['genic'] = all_STRs['downstream_gene_dist'] == 0
for key in 'exonic', 'genic', 'UTR', 'UTR3', 'UTR5', 'promoter':
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

gwas_sig_STR_dfs = []
for phenotype in phenotypes:
    print(f"Loading gwas sig STRs for pheno {phenotype} ... ", flush=True, end='')
    df = pd.read_csv(
        f'{ukb}/association/results/{phenotype}/my_str/results.tab',
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        usecols=('chrom', 'pos', f'p_{phenotype}')
    )
    df = df.loc[df[f'p_{phenotype}'] <= 5e-8, :]
    df = df[['chrom', 'pos']]
    df['gwas_sig'] = True
    gwas_sig_STR_dfs.append(df)
    print('done', flush=True)
gwas_sig_STRs = pd.concat(gwas_sig_STR_dfs).drop_duplicates(['chrom', 'pos']).reset_index(drop=True)
all_STRs = all_STRs.merge(
    gwas_sig_STRs,
    left_on = ['chrom', 'SNPSTR_start_pos'],
    right_on = ['chrom', 'pos'],
    how = 'left',
    suffixes=['', '_other']
)
all_STRs['gwas_sig'] = ~all_STRs['gwas_sig'].isnull()

finemap_STR_dfs = []
for phenotype in phenotypes:
    print(f"Loading FINEMAPed STRs for pheno {phenotype} ... ", flush=True, end='')
    df = pd.read_csv(
        f'{ukb}/finemapping/finemap_results/{phenotype}/summary/all_STR_contribs.tab',
        skiprows = 1,
        delimiter = '\t',
        usecols = ['signal', 'STR', 'pcausal']
    )
    df.rename(columns={'pcausal': 'FINEMAP_pcausal'}, inplace=True)
    df['chrom'] = df['signal'].str.split('_', n=1, expand=True).iloc[:, 0].astype(int)
    df['pos'] = df['STR'].str.split('_', n=1, expand=True).iloc[:, 1].astype(int)
    finemap_STR_dfs.append(df[['chrom', 'pos', 'FINEMAP_pcausal']])
    print('done', flush=True)

finemap_STRs = pd.concat(finemap_STR_dfs).groupby(['chrom', 'pos']).max().reset_index()
all_STRs = all_STRs.merge(
    finemap_STRs,
    left_on = ['chrom', 'SNPSTR_start_pos'],
    right_on = ['chrom', 'pos'],
    how='left',
    suffixes=['', '_other']
)
all_STRs['FINEMAPed'] = False

all_STRs.to_csv(
    f'{annotation_dir}/all_STRs.tab',
    sep='\t',
    index=False
)

# assert only nulls are MHC
'''
assert np.all(all_STRs.loc[all_STRs['gwas_sig'] & all_STRs['FINEMAP_pcausal'].isnull(), 'chrom'] == 6)
assert np.all(all_STRs.loc[all_STRs['gwas_sig'] & all_STRs['FINEMAP_pcausal'].isnull(), 'pos'] >= 25e6)
assert np.all(all_STRs.loc[all_STRs['gwas_sig'] & all_STRs['FINEMAP_pcausal'].isnull(), 'pos'] <= 33.5e6)
'''

# exclude not gwas sig and non MHC
compare_STRs = all_STRs.loc[all_STRs['gwas_sig'] & ~all_STRs['FINEMAP_pcausal'].isnull(), :].copy()
compare_STRs['FINEMAPed'] = compare_STRs['FINEMAP_pcausal'] >= 0.9

# graph number of STRs by FINEMAP p causal
xs = np.arange(0.00, 1.01, 0.01)
fig = bokeh.plotting.figure(
    title=f'Number of sig_gwas STRs by FINEMAP posterior causality',
    y_axis_label = 'Number STRs',
    x_axis_label = 'FINEMAP posterior causality at least this',
    width=1200,
    height=800,
    y_axis_type = 'log'
)
ys = []
for x in xs:
    ys.append(np.sum(compare_STRs['FINEMAP_pcausal'] >= x))

fig.background_fill_color = None
fig.border_fill_color = None
fig.grid.grid_line_color = None
fig.toolbar_location = None
fig.title.text_font_size = '30px'
fig.axis.axis_label_text_font_size = '26px'
fig.axis.major_label_text_font_size = '20px'
fig.line(xs, ys, legend_label='Number STRs')
fig.legend.label_text_font_size = '22px'
bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/FINEMAP_prob_counts.png')
bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/FINEMAP_prob_counts.svg')

with open(f'{ukb}/post_finemapping/results/enrichments.tab', 'w') as out:
    out.write(
        f'From among phenotypes: {phenotypes}\n'
        'All categories are among protein coding genes, except where specified\n'
        'Transcript support level >= 2 required\n'
        'UTR5 and UTR3 must be explicitly labeled as such, not just UTR\n'
        'intronic is not exonic, UTR3, UTR5, but also not coding or UTR (which remove a few more)\n'
        'intergenic is not in a protein coding gene\n'
        'transcribed_non_protein is transcribed and not in a protein coding gene\n'
        'gwas sig is p<=5e-8\n'
        'pcausal cutoff for FINEMAPed is >= 0.9, also required to be gwas sig\n'
    )
    out.write(
        'category\tall_gwas_sig_STRs_count\t#_in_cat\t%_in_cat\tFINEMAPed_STRs_count\t#_in_cat\t%_subset_in_cat\tp_val\n'
    )
    for category in 'exonic', 'UTR5', 'UTR3', 'intronic', 'intergenic', 'transcribed_non_protein', 'eSTR', 'FM_eSTR':
        binary_test(out, compare_STRs, compare_STRs[category], category)
        graph_relative_rate_by_pcausal(compare_STRs, compare_STRs[category], category)
    for period in range(1,7):
        binary_test(out, compare_STRs, compare_STRs['period'] == period, f'period_is_{period}')
        graph_relative_rate_by_pcausal(compare_STRs, compare_STRs['period'] == period, f'period_is_{period}')
    for structure in 'HAIRP', 'QUAD', 'IMOT':
        binary_test(out, compare_STRs, compare_STRs['structure'] == structure, f'structure_is_{structure}')
        graph_relative_rate_by_pcausal(compare_STRs, compare_STRs['structure'] == structure, f'structure_is_{structure}')
    binary_test(out, compare_STRs, (compare_STRs['structure'] != 'UNF') & ~compare_STRs['structure'].isnull(), 'any_stalling_struture')
    graph_relative_rate_by_pcausal(compare_STRs, (compare_STRs['structure'] != 'UNF') & ~compare_STRs['structure'].isnull(), 'any_stalling_struture')

    graph_relative_rate_by_pcausal(compare_STRs, compare_STRs['promoter'], 'promoter')
    repeat_units = sorted(
        set(compare_STRs['canonical_unit']),
        key = lambda seq : (len(seq), seq)
    )
    repeat_units.remove('None')
    for repeat_unit in repeat_units:
        binary_test(out, compare_STRs, compare_STRs['canonical_unit'] == repeat_unit, f'canonical_repeat_unit_is_{repeat_unit}')
    # every STR with more than 1k appearances in all STRs
    for repeat_unit in ('A', 'C', 'AC', 'AT', 'AG', 'AAT', 'AAC', 'AAAT', 'AAAC', 'AGAT', 'AAAG', 'AAGG', 'AATG'):
        graph_relative_rate_by_pcausal(compare_STRs, compare_STRs['canonical_unit'] == repeat_unit, f'unit_is_{repeat_unit}')

    out.write('\n')
    fig, ax = plt.subplots()
    ax.set_xlim(0, np.max(compare_STRs['mean_len']))
    ax.set_title('CDF of STR mean len')
    ax.set_ylabel('density')
    ax.set_xlabel('distance (bp)')
    _BetterCDF(compare_STRs['mean_len'], ax)
    p_val = scipy.stats.mannwhitneyu(compare_STRs['mean_len'], compare_STRs.loc[compare_STRs['FINEMAPed'], 'mean_len'])[1]
    out.write(f'Mann-Whitney p_val mean len FINEMAPed vs all gwas sig: {p_val}\n')
    _BetterCDF(compare_STRs.loc[compare_STRs['FINEMAPed'], 'mean_len'], ax)
    legends = ['gwas sig STRs', 'FINEMAPed']
    ax.text(0.5, -0.1, "Comparing FINEMAPed to gwas_sig. gwas_sig is association p<=5e-8. FINEMAPed is that and also FINEMAP posterior causal >= 0.9. Both groups exclude MHC", ha="center", transform=ax.transAxes)
    ax.text(0.5, -1.1, f"Medians: {np.median(compare_STRs['mean_len'])} (gwas_sig) {np.median(compare_STRs.loc[compare_STRs['FINEMAPed'], 'mean_len'])} (FINEMAPed)", ha="center", transform=ax.transAxes)
    plt.legend(legends)
    plt.savefig(f'{ukb}/post_finemapping/results/mean_len_cdf.png')
    plt.savefig(f'{ukb}/post_finemapping/results/mean_len_cdf.pdf')

    out.write('\n')
    out.write('dist to nearest\tstream\tfrom among\tMann-Whitney p_val FINEMAPed vs all gwas sig\n')
    for category, class_ in ('gene', 'intergenic'), ('exon', 'intronic'):
        for stream in 'upstream', 'downstream':
            col = f'{stream}_{category}_dist'
            subbed_data = compare_STRs.loc[compare_STRs[class_], :]
            fig, ax = plt.subplots()
            ax.set_xlim(0, np.max(subbed_data[col]))
            ax.set_title(f'CDF of {col} from within {class_} STRs')
            ax.set_ylabel('density')
            ax.set_xlabel('distance (bp)')
            _BetterCDF(subbed_data[col], ax)
            p_val = scipy.stats.mannwhitneyu(subbed_data[col], subbed_data.loc[subbed_data['FINEMAPed'], col])[1]
            out.write(f'{category}\t{stream}\t{class_} STRs\t{p_val}\n')
            _BetterCDF(subbed_data.loc[subbed_data['FINEMAPed'], col], ax)
            legends = ['gwas sig STRs', 'FINEMAPed']
            ax.text(0.5, -0.1, "Comparing FINEMAPed to gwas_sig. gwas_sig is association p<=5e-8. FINEMAPed is that and also FINEMAP posterior causal >= 0.9", ha="center", transform=ax.transAxes)
            ax.text(0.5, -1.1, f"Medians: {np.median(subbed_data[col])} (gwas_sig) {np.median(subbed_data.loc[subbed_data['FINEMAPed'], col])} (FINEMAPed)", ha="center", transform=ax.transAxes)
            plt.legend(legends)
            plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf.png')
            plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf.pdf')
            ax.set_xlim(0, 50000)
            subbed_data = subbed_data.loc[subbed_data[col] < 50000, :]
            p_val = scipy.stats.mannwhitneyu(subbed_data[col], subbed_data.loc[subbed_data['FINEMAPed'], col])[1]
            out.write(f'{category}\t{stream}\t{class_} STRs where dist < 50000\t{p_val}\n')
            plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf_50kbp.png')
            plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf_50kbp.pdf')
'''
if __name__ == '__main__':
    main()
'''
