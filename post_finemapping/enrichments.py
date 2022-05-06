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
import polars as pl
import scipy.stats
import scipy.stats.contingency
import seaborn

import annotation_utils
import phenotypes

other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('--load', action='store_true', default=False)
parser.add_argument('--calc', action='store_true', default=False)
args = parser.parse_args()

if args.load:
    # pos (start), snpstr_pos (hipstr)
    all_STRs = pl.read_csv(
        f'{ukb}/snpstr/flank_trimmed_vcf/vars.tab',
        sep='\t'
    )
    # pos (hisptr)
    snpstr_strs = pl.scan_csv(
        f'{ukb}/snpstr/str_loci.txt',
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
    ).drop('snpstr_pos').distinct(
        subset=['chrom', 'pos']
    ).collect()
    assert ~np.any(np.isnan(all_STRs['chrom'].to_numpy()))
    assert ~np.any(np.isnan(all_STRs['pos'].to_numpy()))
    assert ~np.any(np.isnan(all_STRs['end_pos'].to_numpy()))
    assert ~np.any(np.isnan(all_STRs['SNPSTR_start_pos'].to_numpy()))
    # pos (start), end_pos, SNPSTR_start_pos (hipstr)

    repeat_units = pl.read_csv(
        f'{ukb}/snpstr/repeat_units.tab',
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

    '''
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
    '''

    loci_summary_dfs = []
    for chrom in range(1, 23):
        distribution_stats = pl.read_csv(
            f'{ukb}/export_scripts/intermediate_results/chr{chrom}_loci_summary.tab',
            sep='\t',
        )
        loci_summary_dfs.append(distribution_stats)
    loci_summaries = pl.concat(loci_summary_dfs)
    n_before = all_STRs.shape[0]
    all_STRs = all_STRs.join(
        loci_summaries,
        how='left',
        left_on=['chrom', 'SNPSTR_start_pos'],
        right_on=['chr', 'pos'],
        suffix='_other'
    )
    assert n_before == all_STRs.shape[0]
    print('Calculating mean lens ... ', flush=True, end='')
    all_STRs = all_STRs.with_column(
        pl.Series([sum(key*val for (key, val) in ast.literal_eval(allele_dist).items()) for allele_dist in all_STRs['allele_dist']]).alias('mean_len')
    )
    print('done', flush=True)

    print('Loading eSTRs ... ')
    eSTRs = pl.read_csv(
        f'{ukb}/misc_data/eSTR/eSTRs.csv',
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
        f'{ukb}/misc_data/gencode/gencode.v38lift37.annotation.without_chr.sorted.gene.gff3',
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
    all_STRs['promoter'] = False
    for row in range(genes.shape[0]):
        if not 'gene_type=protein_coding' in genes[row, 'kvps']:
            continue
        if genes[row, 'chrom'] not in set(str(x) for x in range(1, 23)):
            continue
        tss = genes[row, 'TSS']
        all_STRs.loc[
            (all_STRs['chrom'] == int(genes[row, 'chrom'])) & (
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
        f'{ukb}/post_finemapping/intermediate_results/concordant_causal_STR_candidates.tab',
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
    #for phenotype in ['albumin', 'cystatin_c']: #phenotypes.phenotypes_in_use:
    for phenotype in phenotypes.phenotypes_in_use:
        print(phenotype, flush=True)
        gwsig_df = pl.scan_csv(
            f'{ukb}/association/results/{phenotype}/my_str/results.tab',
            sep = '\t',
        ).filter(pl.col(f'p_{phenotype}') < 5e-8).select(['chrom', 'pos']).collect()
        gwsig_loci = gwsig_loci.union(zip(gwsig_df['chrom'].to_numpy(), gwsig_df['pos'].to_numpy()))
    gwsig_loci = [f'{chrom}_{pos}' for (chrom, pos) in gwsig_loci]
    print('matching gwsig loci .... ', flush=True)
    all_STRs = all_STRs.with_column(
        (pl.col('chrom').cast(str) + '_' +  pl.col('SNPSTR_start_pos').cast(str)).is_in(gwsig_loci).alias('gwsig')
    )

    '''
    finemap_STRs = pd.concat(finemap_STR_dfs).groupby(['chrom', 'pos']).max().reset_index()
    all_STRs = all_STRs.merge(
        finemap_STRs,
        left_on = ['chrom', 'SNPSTR_start_pos'],
        right_on = ['chrom', 'pos'],
        how='left',
        suffixes=['', '_other']
    )
    '''

    '''
    assert np.all(all_STRs.loc[all_STRs['gwas_sig'] & all_STRs['FINEMAP_pcausal'].isnull(), 'chrom'] == 6)
    assert np.all(all_STRs.loc[all_STRs['gwas_sig'] & all_STRs['FINEMAP_pcausal'].isnull(), 'pos'] >= 25e6)
    assert np.all(all_STRs.loc[all_STRs['gwas_sig'] & all_STRs['FINEMAP_pcausal'].isnull(), 'pos'] <= 33.5e6)
    '''

    all_STRs.write_csv(
        f'{ukb}/post_finemapping/intermediate_results/enrichments_df.tab',
        sep='\t'
    )
    exit()

if args.calc:
    print('Reading CSV ...', flush=True)
    compare_STRs = pl.scan_csv(
        f'{ukb}/post_finemapping/intermediate_results/enrichments_df.tab',
        sep='\t',
        dtypes={'is_causal_STR_candidate': bool}
    ).with_column(
        (~pl.col('is_causal_STR_candidate').is_null()).alias('is_causal_STR_candidate')
    ).collect()

    # assert only nulls are MHC

    def compare_categories(out, compare_STRs, category, category_string):
        print(f'Comparing {category_string} ...', flush=True)
        n_strs = compare_STRs.shape[0]
        n_cat = compare_STRs.filter(category).shape[0]
        n_strs_gwsig = compare_STRs.filter('gwsig').shape[0]
        n_cat_gwsig = compare_STRs.filter(category & pl.col('gwsig')).shape[0]
        #n_strs_gwsig = compare_STRs.filter(pl.col('p_val') <= 5e-8).shape[0]
        #n_cat_gwsig = compare_STRs.filter(category & (pl.col('p_val') <= 5e-8)).shape[0]
        assert n_strs - n_cat - n_strs_gwsig + n_cat_gwsig == compare_STRs.filter(~category & ~pl.col('gwsig')).shape[0]
        n_strs_fm= compare_STRs.filter('is_causal_STR_candidate').shape[0]
        n_cat_fm = compare_STRs.filter(category & pl.col('is_causal_STR_candidate')).shape[0]
        p_vals = []
        for table in (
            [[n_strs - n_cat - n_strs_gwsig + n_cat_gwsig, n_strs_gwsig - n_cat_gwsig], [n_cat - n_cat_gwsig, n_cat_gwsig]],
            [[n_strs - n_cat - n_strs_fm + n_cat_fm, n_strs_fm - n_cat_fm], [n_cat - n_cat_fm, n_cat_fm]],
            [[n_strs_gwsig - n_cat_gwsig - n_strs_fm + n_cat_fm, n_strs_fm - n_cat_fm], [n_cat_gwsig - n_cat_fm, n_cat_fm]],
        ):
            if np.any(np.array(table) < 5):
                p_vals.append(scipy.stats.fisher_exact(table)[1])
            else:
                p_vals.append(scipy.stats.chi2_contingency(table, correction=False)[1])

        out.write(
            f'{category_string}\t'
            f'{n_cat}/{n_strs} ({n_cat/n_strs*100:.4f}%)\t'
            f'{n_cat_gwsig}/{n_strs_gwsig} ({n_cat_gwsig/n_strs_gwsig*100:.4f}%)\t'
            f'{n_cat_fm}/{n_strs_fm} ({n_cat_fm/n_strs_fm*100:.4f}%)\t'
            f'{p_vals[0]}\t'
            f'{p_vals[1]}\t'
            f'{p_vals[2]}\n'
        )
        out.flush()

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
    '''

    print('Running calcs ... ', flush=True)
    with open(f'{ukb}/post_finemapping/results/enrichments.tab', 'w') as out:
        out.write(
            'category\t'
            'n_cat/n_strs (%)\t'
            'n_cat_gwsig/n_strs_gwsig (%)\t'
            'n_cat_fm/n_strs_fm (%)\t'
            'p_all_strs_vs_gwsig\t'
            'p_all_strs_vs_causal\t'
            'p_gwsig_vs_causal\n'
        )
        out.flush()
        compare_categories(out, compare_STRs, pl.col('exonic') & (pl.col('period') == 3), 'exonic_trinucs')
        compare_categories(out, compare_STRs, ~pl.col('exonic') & (pl.col('period') == 3), 'non_exonic_trinucs')
        compare_categories(out, compare_STRs, pl.col('exonic') & (pl.col('period') != 3), 'exonic_not_trinucs')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'A', 'polyA')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'C', 'polyC')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AC', 'polyAC')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AT', 'polyAT')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AG', 'polyAG')
        compare_categories(out, compare_STRs, ~pl.col('canonical_unit').is_in(['AC', 'AT', 'AG']) & (pl.col('period') == 2), 'other_dinucs')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'CCG', 'polyCCG')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAT', 'polyAAAT')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAC', 'polyAAAC')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AGAT', 'polyAGAT')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAG', 'polyAAAG')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAGG', 'polyAAGG')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AATG', 'polyAATG')
        compare_categories(out, compare_STRs, ~pl.col('canonical_unit').is_in(['AAAT', 'AAAC', 'AGAT',' AAAG', 'AAGG', 'AATG']) & (pl.col('period') == 4), 'other_tetranucs')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAAC', 'polyAAAAC')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAAT', 'polyAAAAT')
        compare_categories(out, compare_STRs, ~pl.col('canonical_unit').is_in(['AAAAC', 'AAAAT']) & (pl.col('period') == 5), 'other_pentanucs')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAAAC', 'polyAAAACC')
        compare_categories(out, compare_STRs, ~pl.col('canonical_unit').is_in(['AAAAAC']) & (pl.col('period') == 6), 'other_hexanucs')
        compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'None', 'unclear_repeat_unit')
        compare_categories(out, compare_STRs, pl.col('coding'), 'coding')
        compare_categories(out, compare_STRs, pl.col('UTR3'), "3'-UTR")
        compare_categories(out, compare_STRs, pl.col('UTR5'), "5'-UTR")
        compare_categories(out, compare_STRs, pl.col('intronic'), 'intronic')
        compare_categories(out, compare_STRs, pl.col('genic'), 'genic')
        compare_categories(out, compare_STRs, pl.col('transcribed') & ~pl.col('genic'), 'transcribed_non_protein')
        compare_categories(out, compare_STRs, pl.col('promoter') & ~pl.col('genic'), 'upstream_promoter')
        compare_categories(out, compare_STRs, pl.col('eSTR'), 'eSTR')
        compare_categories(out, compare_STRs, pl.col('FM_eSTR'), 'FM_eSTR')

    exit()

def barplot_fig(cats, data, ps):
    fig = bokeh.plotting.figure(
        y_axis_label = 'Prevalence (%)',
        x_range = bokeh.models.FactorRange(*cats),
        width=len(data)//3*200,
        height=800,
        toolbar_location=None
    )
    arr_data = np.array(data)
    arr_data[arr_data == 0] = np.nan
    fig.vbar(x=cats, top=arr_data, width=0.9)
    fig.background_fill_color = None
    fig.border_fill_color = None
    fig.x_range.range_padding = 0.1
    fig.xaxis.major_label_orientation = 1
    fig.xgrid.grid_line_color = None
    y_range = max(data) # since plots are relative to zero
    for group in range(len(cats)//3):
        compare_bars(fig, y_range, cats[group*3], cats[group*3+1], data[group*3], data[group*3+1], data[group*3:group*3+2], ps[group*3], 0, -1, 1)
        compare_bars(fig, y_range, cats[group*3+1], cats[group*3+2], data[group*3+1], data[group*3+2], data[group*3+1:group*3+3], ps[group*3+2], 1, 0, 1)
        compare_bars(fig, y_range, cats[group*3], cats[group*3+2], data[group*3], data[group*3+2], data[group*3:group*3+3], ps[group*3+1], -1, 1, 2)
    return fig

def compare_bars(fig, y_range, x1, x2, y1, y2, ys, p_val, x1step, x2step, ystep):
    x1 = [*x1, .15*x1step]
    x2 = [*x2, .15*x2step]
    y_step_size = 0.05*y_range
    top = max(ys) + ystep*y_step_size
    fig.line(x=[x1, x1], y=[y1, top])
    fig.line(x=[x2, x2], y=[y2, top])
    fig.line(x=[x1, x2], y=[top, top])
    text=f'(p={p_val:.1g})'
    if p_val < 0.05/50:
        text = '** ' + text
    elif p_val < 0.05:
        text = '* ' + text

    cds = bokeh.models.ColumnDataSource(dict(
        x=[x1[:2]],
        y=[top+0.01*y_step_size],
        text=[text]
    ))
    fig.add_layout(bokeh.models.LabelSet(x='x', y='y', x_offset=.05, text='text', source=cds, text_font_size='12px'))

datas = {}
for fig_loc in 'main', 'supp':
    datas[fig_loc] = {}
    for bar_loc in 'left', 'right':
        datas[fig_loc][bar_loc] = {}
        datas[fig_loc][bar_loc]['data'] = []
        datas[fig_loc][bar_loc]['cat'] = []
        datas[fig_loc][bar_loc]['ps'] = []

with open(f'{ukb}/post_finemapping/results/enrichments.tab') as stats:
    next(stats) #skip header
    for line in stats:
        if 'CCG' in line:
            continue
        split = line.split()
        category = split[0]
        ori_perc = float(split[2][1:-2])
        if 'other' in category or 'poly' in category or 'unclear' in category:
            fig_loc = 'supp'
        else:
            fig_loc = 'main'
        if ori_perc >= 10:
            bar_loc = 'right'
        else:
            bar_loc = 'left'
        curr_data = datas[fig_loc][bar_loc]['data']
        curr_cat = datas[fig_loc][bar_loc]['cat']
        curr_ps = datas[fig_loc][bar_loc]['ps']
        curr_cat.extend([
            (category, f'all STRs'),
            (category, f'gwsig STRs'),
            (category, 'high-confidence STRs')
        ])
        curr_ps.extend([float(split[idx]) for idx in range(7, 10)])

        curr_data.append(ori_perc)
        curr_data.append(float(split[4][1:-2]))
        curr_data.append(float(split[6][1:-2]))

fig = bokeh.layouts.row([
    barplot_fig(datas['main']['left']['cat'], datas['main']['left']['data'], datas['main']['left']['ps']),
    barplot_fig(datas['main']['right']['cat'], datas['main']['right']['data'], datas['main']['right']['ps']),
])

bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/enrichment_barplots.svg')
bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/enrichment_barplots.png')

fig = bokeh.layouts.row([
    barplot_fig(datas['supp']['left']['cat'], datas['supp']['left']['data'], datas['supp']['left']['ps']),
    barplot_fig(datas['supp']['right']['cat'], datas['supp']['right']['data'], datas['supp']['right']['ps']),
])
bokeh.io.export_svg(fig, filename=f'{ukb}/post_finemapping/results/enrichment_barplots_supp.svg')
bokeh.io.export_png(fig, filename=f'{ukb}/post_finemapping/results/enrichment_barplots_supp.png')
exit()

# old code for continuous comparisons
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

with open(f'{ukb}/post_finemapping/results/enrichments.tab', 'w') as out:
    f, axs = plt.subplots(1,2, figsize=(9,5), sharey=True)
    seaborn.boxplot(
        y='mean_len',
        data=compare_STRs.filter('concordantly_finemapped').select('mean_len').to_pandas(),
        ax=axs[0],
        fliersize=0,
        boxprops={'facecolor': 'none'},
        color=None
    )
    seaborn.swarmplot(
        y='mean_len',
        data=compare_STRs.filter('concordantly_finemapped').select('mean_len').to_pandas(),
        ax=axs[0],
        color='black'
    ).set(title='concordantly finemapped STRs')
    ax = seaborn.violinplot(
        y='mean_len',
        data=compare_STRs.select('mean_len').to_pandas(),
        inner=None,
        ax=axs[1]
    )
    ax.set(title='all STRs in finemapping regions with p-value <= 1e-10')
    ax.collections[0].set_facecolor((0,0,0,0))
    seaborn.boxplot(
        y='mean_len',
        data=compare_STRs.select('mean_len').to_pandas(),
        ax=axs[1],
        boxprops={'facecolor': 'none'},
        fliersize=0
    )
    plt.savefig(f'{ukb}/post_finemapping/results/mean_len_comparison.png')
    p_val = scipy.stats.mannwhitneyu(compare_STRs['mean_len'].to_numpy(), compare_STRs.filter('concordantly_finemapped')['mean_len'].to_numpy())[1]
    out.write(f'Mann-Whitney p_val mean len concordantly fine-mapped vs all gwas sig in finemapping regions: {p_val}\n')

    f, axs = plt.subplots(1,2, figsize=(9,5), sharey=True)
    seaborn.boxplot(
        y='upstream_gene_dist',
        data=compare_STRs.filter(pl.col('concordantly_finemapped') & pl.col('intergenic') & (pl.col('upstream_gene_dist') <= 50_000)).select('upstream_gene_dist').to_pandas(),
        ax=axs[0],
        fliersize=0,
        boxprops={'facecolor': 'none'},
        color=None
    )
    seaborn.swarmplot(
        y='upstream_gene_dist',
        data=compare_STRs.filter(pl.col('concordantly_finemapped') & pl.col('intergenic') & (pl.col('upstream_gene_dist') <= 50_000)).select('upstream_gene_dist').to_pandas(),
        ax=axs[0],
        color='black'
    ).set(title='concordantly finemapped STRs')
    ax = seaborn.violinplot(
        y='upstream_gene_dist',
        data=compare_STRs.filter(pl.col('intergenic') & (pl.col('upstream_gene_dist') <= 50_000)).select('upstream_gene_dist').to_pandas(),
        inner=None,
        ax=axs[1]
    )
    ax.set(title='all STRs in finemapping regions with p-value <= 1e-10')
    ax.collections[0].set_facecolor((0,0,0,0))
    seaborn.boxplot(
        y='upstream_gene_dist',
        data=compare_STRs.filter(pl.col('intergenic') & (pl.col('upstream_gene_dist') <= 50_000)).select('upstream_gene_dist').to_pandas(),
        ax=axs[1],
        boxprops={'facecolor': 'none'},
        fliersize=0
    )
    plt.savefig(f'{ukb}/post_finemapping/results/intergenic_upstream_gene_dist_comparison.png')
    p_val = scipy.stats.mannwhitneyu(
        compare_STRs.filter(pl.col('intergenic') & (pl.col('upstream_gene_dist') <= 50_000))['upstream_gene_dist'].to_numpy(),
        compare_STRs.filter(pl.col('concordantly_finemapped') & pl.col('intergenic') & (pl.col('upstream_gene_dist') <= 50_000))['upstream_gene_dist'].to_numpy(),
    )[1]
    out.write(f'Mann-Whitney p_val dist to nearest upstream gene in intergenic STRs (capped at 50kbp), concordantly fine-mapped vs all gwas sig in finemapping regions: {p_val}\n')


    out.write('\n')
    out.write('dist to nearest\tstream\tfrom among\tMann-Whitney p_val FINEMAPed vs all gwas sig\n')
    for category, class_ in [('gene', 'intergenic'), ('exon', 'intronic')]:
        for stream in ['upstream', 'downstream']:
            col = f'{stream}_{category}_dist'
            subbed_data = compare_STRs.filter(pl.col(class_))
            fig, ax = plt.subplots()
            ax.set_xlim(0, np.max(subbed_data[col].to_numpy()))
            ax.set_title(f'CDF of {col} from within {class_} STRs')
            ax.set_ylabel('density')
            ax.set_xlabel('distance (bp)')
            _BetterCDF(subbed_data[col].to_numpy(), ax)
            '''
            p_val = scipy.stats.mannwhitneyu(subbed_data[col], subbed_data.loc[subbed_data['FINEMAPed'], col])[1]
            out.write(f'{category}\t{stream}\t{class_} STRs\t{p_val}\n')
            '''
            _BetterCDF(subbed_data.filter(pl.col('concordantly_finemapped'))[col].to_numpy(), ax)
            legends = ['gwas sig STRs', 'concordantly_finemapped']
            '''
            ax.text(0.5, -0.1, "Comparing FINEMAPed to gwas_sig. gwas_sig is association p<=5e-8. FINEMAPed is that and also FINEMAP posterior causal >= 0.9", ha="center", transform=ax.transAxes)
            ax.text(0.5, -1.1, f"Medians: {np.median(subbed_data[col])} (gwas_sig) {np.median(subbed_data.loc[subbed_data['FINEMAPed'], col])} (FINEMAPed)", ha="center", transform=ax.transAxes)
            '''
            plt.legend(legends)
            plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf.png')
            plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf.pdf')
            ax.set_xlim(0, 50000)
            '''
            subbed_data = subbed_data.loc[subbed_data[col] < 50000, :]
            p_val = scipy.stats.mannwhitneyu(subbed_data[col], subbed_data.loc[subbed_data['FINEMAPed'], col])[1]
            out.write(f'{category}\t{stream}\t{class_} STRs where dist < 50000\t{p_val}\n')
            '''
            plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf_50kbp.png')
            plt.savefig(f'{ukb}/post_finemapping/results/{col}_cdf_50kbp.pdf')
