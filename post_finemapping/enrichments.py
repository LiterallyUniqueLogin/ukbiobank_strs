#!/usr/bin/env python3

import argparse
import os
from typing import List

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats

import annotation_utils

ukb = os.environ['UKB']

'''
ALL STRs
And with SIG <= 5e-8
And with FINEMAP pcausal >= .9
'''
def main():
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
    all_STRs = all_STRs[['chrom', 'pos', 'snpstr_pos']]
    all_STRs.rename(columns={'snpstr_pos': 'SNPSTR_start_pos'}, inplace=True)
    all_STRs = all_STRs.drop_duplicates(subset=['chrom', 'pos']).reset_index(drop=True)
    # pos (start), SNPSTR_start_pos (hipstr)

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
    for key in 'exonic', 'genic', 'UTR', 'UTR3', 'UTR5':
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

    all_STRs.to_csv(
        f'{annotation_dir}/all_STRs.tab',
        sep='\t',
        index=False
    )

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
        how = 'left'
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
        df = df.loc[df['pcausal'] >= 0.9, :]

        df['chrom'] = df['signal'].str.split('_', n=1, expand=True).iloc[:, 0].astype(int)
        df['pos'] = df['STR'].str.split('_', n=1, expand=True).iloc[:, 1].astype(int)
        df['FINEMAPed'] = True
        finemap_STR_dfs.append(df[['chrom', 'pos', 'FINEMAPed']])
        print('done', flush=True)

    finemap_STRs = pd.concat(finemap_STR_dfs).drop_duplicates(['chrom', 'pos']).reset_index(drop=True)
    all_STRs = all_STRs.merge(
        finemap_STRs,
        left_on = ['chrom', 'SNPSTR_start_pos'],
        right_on = ['chrom', 'pos'],
        how='left'
    )
    all_STRs['FINEMAPed'] = ~all_STRs['FINEMAPed'].isnull()

    n_strs = all_STRs.shape[0]
    with open(f'{ukb}/post_finemapping/results/binary_enrichments.tab', 'w') as out:
        out.write(
            f'From among phenotypes: {phenotypes}\n'
            'All categories are among protein coding genes, except where specified\n'
            'Transcript support level >= 2 required\n'
            'UTR5 and UTR3 must be explicitly labeled as such, not just UTR\n'
            'intronic is not exonic, UTR3, UTR5, but also not coding or UTR (which remove a few more)\n'
            'intergenic is not in a protein coding gene\n'
            'transcribed_non_protein is transcribed and not in a protein coding gene\n'
        )
        out.write(
            'category\tsubset\tall_STRs_count\t%_in_cat\tsubset_STRs_count\t%_subset_in_cat\tp_val\n'
        )
        for category in 'exonic', 'UTR5', 'UTR3', 'intronic', 'intergenic', 'transcribed_non_protein':
            n_cat = np.sum(all_STRs[category])
            for subset in 'gwas_sig', 'FINEMAPed':
                n_neither = np.sum(~all_STRs[subset] & ~all_STRs[category])
                n_subset = np.sum(all_STRs[subset])
                n_subset_not_cat = np.sum(all_STRs[subset] & ~all_STRs[category])
                n_cat_not_subset = np.sum(~all_STRs[subset] & all_STRs[category])
                n_subset_cat = np.sum(all_STRs[subset] & all_STRs[category])
                p_val = scipy.stats.chi2_contingency([[n_neither, n_subset_not_cat], [n_cat_not_subset, n_subset_cat]])[1]
                out.write(f'{category}\t{subset}\t{n_strs}\t{n_cat/n_strs*100:.4f}%\t{n_subset}\t{n_subset_cat/n_subset*100:.4f}%\t{p_val}\n')

        out.write('\n')
        out.write('dist to nearest\tstream\tfrom among\tcompared to\tp_val\n')
        for category, class_ in ('gene', 'intergenic'), ('exon', 'intronic'):
            for stream in 'upstream', 'downstream':
                col = f'{stream}_{category}_dist'
                subbed_data = all_STRs.loc[all_STRs[class_], :]
                fig, ax = plt.subplots()
                ax.set_xlim(0, np.max(subbed_data[col]))
                ax.set_title(f'CDF of {col} from within {class_} STRs')
                ax.set_ylabel('density')
                ax.set_xlabel('distance (bp)')
                legends = ['all_STRs']
                _BetterCDF(subbed_data[col], ax)
                for subset in 'gwas_sig', 'FINEMAPed':
                    p_val = scipy.stats.mannwhitneyu(subbed_data[col], subbed_data.loc[subbed_data[subset], col])[1]
                    out.write(f'{category}\t{stream}\t{class_} STRs\t{subset} STRs\t{p_val}\n')
                    _BetterCDF(subbed_data.loc[subbed_data[subset], col], ax)
                    legends.append(subset)
                plt.legend(legends)
                plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf.png')
                plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf.svg')

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

if __name__ == '__main__':
    main()
