#!/usr/bin/env python3

import argparse

import numpy as np
import polars as pl
import scipy.stats
import scipy.stats.contingency

other_ethnicities = ['black', 'south_asian', 'chinese', 'irish', 'white_other']

parser = argparse.ArgumentParser()
parser.add_argument('outdir')
parser.add_argument('enrichments_df')
args = parser.parse_args()

print('Reading CSV ...', flush=True)
compare_STRs = pl.scan_csv(
    #f'{ukb}/post_finemapping/intermediate_results/enrichments_df.tab',
    args.enrichments_df,
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
    #n_strs_gwsig = compare_STRs.filter(pl.col('p_val') < 5e-8).shape[0]
    #n_cat_gwsig = compare_STRs.filter(category & (pl.col('p_val') < 5e-8)).shape[0]
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

print('Running calcs ... ', flush=True)
#with open(f'{ukb}/post_finemapping/results/enrichments.tab', 'w') as out:
with open(f'{args.outdir}/enrichments.tab', 'w') as out:
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
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'A', 'poly-A')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'C', 'poly-C')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AC', 'poly-AC')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AG', 'poly-AG')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AT', 'poly-AT')
    compare_categories(out, compare_STRs, ~pl.col('canonical_unit').is_in(['AC', 'AT', 'AG']) & (pl.col('period') == 2), 'other_dinucs')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'CCG', 'poly-CCG')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAT', 'poly-AAAT')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAC', 'poly-AAAC')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AGAT', 'poly-AGAT')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAG', 'poly-AAAG')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAGG', 'poly-AAGG')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AATG', 'poly-AATG')
    compare_categories(out, compare_STRs, ~pl.col('canonical_unit').is_in(['AAAT', 'AAAC', 'AGAT',' AAAG', 'AAGG', 'AATG']) & (pl.col('period') == 4), 'other_tetranucs')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAAC', 'poly-AAAAC')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAAT', 'poly-AAAAT')
    compare_categories(out, compare_STRs, ~pl.col('canonical_unit').is_in(['AAAAC', 'AAAAT']) & (pl.col('period') == 5), 'other_pentanucs')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'AAAAAC', 'poly-AAAACC')
    compare_categories(out, compare_STRs, ~pl.col('canonical_unit').is_in(['AAAAAC']) & (pl.col('period') == 6), 'other_hexanucs')
    compare_categories(out, compare_STRs, pl.col('canonical_unit') == 'None', 'unclear_repeat_unit')
    compare_categories(out, compare_STRs, pl.col('coding'), 'coding')
    compare_categories(out, compare_STRs, pl.col('UTR3'), "3'-UTR")
    compare_categories(out, compare_STRs, pl.col('UTR5'), "5'-UTR")
    compare_categories(out, compare_STRs, pl.col('intronic'), 'intronic')
    compare_categories(out, compare_STRs, pl.col('genic'), 'genic')
    compare_categories(out, compare_STRs, pl.col('transcribed') & ~pl.col('genic'), 'transcribed_non_protein')
    compare_categories(out, compare_STRs, pl.col('upstream_promoter'), 'upstream_promoter')
    compare_categories(out, compare_STRs, pl.col('eSTR'), 'eSTR')
    compare_categories(out, compare_STRs, pl.col('FM_eSTR'), 'FM_eSTR')

