#!/usr/bin/env python3

import os

import matplotlib.cm
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from qqman import qqman

ukb = os.environ['UKB']

def validate_our_code():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))

    print('loading csvs ... ', flush=True)
    my_df = pd.read_csv(
        f'{ukb}/association/plots/input/eosinophil_count/my_imputed_snp_chr21_results.tab',
        header=0,
        delimiter='\t',
        usecols=['chrom', 'pos', 'alleles', 'p_eosinophil_count', 'locus_filtered'],
        dtype={'chrom': int, 'pos': int, 'alleles': str, 'p_eosinophil_count': float, 'locus_filtered': str},
    )
    my_df.rename(columns={'chrom': 'CHR', 'pos': 'BP', 'p_eosinophil_count': 'P'}, inplace=True)
    alleles = my_df['alleles'].str.split(',', 1, expand=True)
    print(my_df.shape)
    
    my_df = my_df[my_df['locus_filtered'] == 'False']
    print(my_df.shape)
    my_df['REF'] = alleles.iloc[:, 0]
    my_df['ALT'] = alleles.iloc[:, 1]

    plink_df = pd.read_csv(
        f'{ukb}/association/results/eosinophil_count/plink_snp/results.tab',
        header=0,
        delimiter='\t',
        usecols=['#CHROM', 'POS', 'P', 'REF', 'ALT'],
        dtype={'#CHROM': int, 'POS': int, 'P': float}
    )
    plink_df.rename(columns={'#CHROM': 'CHR', 'POS': 'BP'}, inplace=True)
    plink_df = plink_df[plink_df['CHR'] == 21]
    print(plink_df.shape)

    plink_df = plink_df.merge(
        my_df,
        how='inner',
        on=['BP', 'REF', 'ALT'],
        suffixes=(None, '_mine')
    )
    plink_df = plink_df[['CHR', 'BP', 'P', 'REF', 'ALT']]
    print(plink_df.shape)

    my_df = my_df.merge(
        plink_df,
        how='inner',
        on=['BP', 'REF', 'ALT'],
        suffixes=(None, '_theirs')
    )
    my_df = my_df[['CHR', 'BP', 'P']]
    print(my_df.shape)

    print('plotting manhattans ... ', flush=True)
    qqman.manhattan(
        my_df,
        cmap=matplotlib.cm.get_cmap('viridis'),
        ax=ax1,
        title="P-values from this study's code",
        suggestiveline=False
    )

    qqman.manhattan(
        plink_df,
        cmap=matplotlib.cm.get_cmap('viridis'),
        ax=ax2,
        title='Plink p-values',
        suggestiveline=False
    )

    plt.savefig(f'{ukb}/export_scripts/results/validate_our_code.png')

# currently unused, doesn't seem like a figure we want
def compare_to_panukbb_b():
    print('loading csv ... ', flush=True)
    df = pd.read_csv(
        f'{ukb}/misc_data/snp_summary_stats/bilirubin/neale/biomarkers-30840-both_sexes-irnt.tsv.bgz',
        compression='gzip',
        header=0,
        delimiter='\t',
        usecols=['chr', 'pos', 'pval_EUR'],
        dtype={'chr': str, 'pos': int, 'pval_EUR': float}
    )
    df.rename(columns={'chr': 'CHR', 'pos': 'BP', 'pval_EUR': 'P'}, inplace=True)
    df = df[df['CHR'] != 'X']
    df = df.astype({'CHR': int})
    df.loc[df['P'] < 1e-100, 'P'] = 1e-100

    print('plotting manhattan ... ', flush=True)
    qqman.manhattan(
        df,
        f'{ukb}/export/figures/compare_to_pan_ukbb_b.png',
        suggestiveline=False
    )

def scatter_with_panukbb():
    print('loading panukbb ... ', flush=True)
    panukbb_df = pd.read_csv(
        f'{ukb}/misc_data/snp_summary_stats/bilirubin/neale/biomarkers-30840-both_sexes-irnt.tsv.bgz',
        compression='gzip',
        header=0,
        delimiter='\t',
        usecols=['chr', 'pos', 'ref', 'alt', 'pval_EUR'],
        dtype={'chr': str, 'pos': int, 'ref': str, 'alt': str, 'pval_EUR': float}
    )
    panukbb_df.rename(columns={'pval_EUR': 'p'}, inplace=True)
    panukbb_df = panukbb_df[panukbb_df['chr'] != 'X']
    panukbb_df = panukbb_df.astype({'chr': int})

    print('loading plink ... ', flush=True)
    my_pipeline_df = pd.read_csv(
        f'{ukb}/association/results/total_bilirubin/plink_snp/results.tab',
        header=0,
        delimiter='\t',
        usecols=['#CHROM', 'POS', 'P', 'REF', 'ALT'],
        dtype={'#CHROM': int, 'POS': int, 'REF': str, 'ALT': str, 'P': float}
    )
    my_pipeline_df.rename(
        columns={'#CHROM': 'chr', 'POS': 'pos', 'REF': 'ref', 'ALT': 'alt', 'P': 'p'},
        inplace=True
    )

    print('merging ... ', flush=True)
    merged = panukbb_df.merge(
        my_pipeline_df,
        how='inner',
        on=['chr', 'pos', 'ref', 'alt'],
        suffixes=('_panukbb', '_mine')
    )

    merged['p_panukbb'] = np.maximum(1e-300, merged['p_panukbb'])
    merged['p_mine'] = np.maximum(1e-300, merged['p_mine'])

    merged['p_panukbb'] = -np.log10(merged['p_panukbb'])
    merged['p_mine'] = -np.log10(merged['p_mine'])

    merged['p_panukbb'] = np.minimum(50, merged['p_panukbb'])
    merged['p_mine'] = np.minimum(50, merged['p_mine'])

    print('plotting ... ', flush=True)
    figure, ax = plt.subplots(1, 1)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('-log10(p panukbb)')
    ax.set_ylabel('-log10(p our pipeline)')
    # From colorbrewer
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(name='meh', colors=['#e5f5f9', '#99d8c9', '#2ca25f'])
    hexbin = ax.hexbin(merged['p_panukbb'], merged['p_mine'], gridsize=100, cmap=cmap, bins='log')
    cb = figure.colorbar(hexbin, ax=ax)
    cb.set_label('counts')
    figure.show()

    plt.savefig(f'{ukb}/export_scripts/results/panukbb_scatter.png')

if __name__ == '__main__':
    validate_our_code()
    scatter_with_panukbb()
