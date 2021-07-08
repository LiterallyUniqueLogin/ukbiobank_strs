#!/usr/bin/env python3

import argparse
import ast
import os

import cyvcf2
import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import WLS
from statsmodels.stats.weightstats import DescrStatsW

import python_array_utils as utils

ukb = os.environ['UKB']

def prep_dict(d):
    return {float(k): v for k, v in d.items() if not np.isnan(v)}

def rename_column_pd(df, old_name, new_name):
    df.rename(columns = {old_name: new_name}, inplace = True)

def main(readme, phenotype, previous_STRs):
    with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt') as unitfile:
        unit = next(unitfile).strip()

    col_dphen_unit = 'Δphenotype_per_additional_repeat_unit'
    col_dphen_sd = 'Δphenotype_per_s.d._increase_in_repeat_size'

    # ordered list of pairs (colname, description)
    columns = [
        ('chrom', None),
        ('signal_region', '{start bp}_{end bp} of region of association the STR resides in, '
         '1-based, inclusive'),
        ('start_pos', '1-based, inclusive'), #
        ('end_pos', '1-based, inclusive'), #
        ('SNPSTR_start_pos', 'the start position of the STR in the SNP-STR reference panel '
         '(may be smaller that start_pos if HipSTR called this STR with a physically phased '
         'SNP upstream of the STR)'),
        ('SNPSTR_ID', 'ID of the STR in the SNPSTR reference panel'),
        ('period', 'as specified in the SNPSTR reference panel'),
        ("repeat_unit", "as inferred by TRTools from the STR sequence and "
         "the period"),
        ('alleles', 'possible # copies of the repeat unit'),
        ('reference_allele', None),
        ('total_per_allele_dosages', 'sum of imputed dosage of each allele across '
         'both chromsomes in all samples'),
        ('total_hardcall_alleles', 'number of each allele in the population using '
         "Beagle's called most probable phased genotypes"),
        ('total_hardcall_genotypes', 'number of each unphased genotype in the '
         "population using Beagle's called most probable phased genotypes"),
        ('subset_total_per_allele_dosages', 'as total_per_allele_dosages, but '
         'calculated on the subset of samples used in the association'),
        ('subset_total_hardcall_alleles', 'as total_hardcall_alleles , but '
         'calculated on the subset of samples used in the association'),
        ('subset_total_hardcall_genotypes', 'as total_hardcall_genotypes, but '
         'calculated on the subset of samples used in the association'),
        ('subset_multiallelicness', '% of dosage attributed to alleles that are not '
         'the first or second most frequent alleles, by dosage, within the subset of '
         'samples used in the association'),
        ('subset_heterozygosity', 'the expected percent of heterozygous samples in the '
         'association sample subset as calculated from allele frequencies'),
        ('subset_entropy', 'the entropy in the association sample subset calculated '
         'from allele frequencies'),
        ('subset_HWEP', 'the Hardy-Weinberg p-value in the association sample subset '
         'comparing hardcall genotypes and allele frequencies'),
        ('subset_allele_dosage_r2', 'a metric of imputation accuracy, the per-allele R^2 '
         'between dosage-weighted allele length and allele hardcalls across samples'),
        ('association_p_value', 'linear regression of rank inverse normalized phenotype '
         "values vs repeat length on QC'ed sample subset, with covariates"),
        (col_dphen_unit, 'linear regression on raw phenotypes vs repeat length on '
         f"QC'ed sample subset, no covariates included, phenotype measured in {unit}"),
        (col_dphen_sd, 'linear regression on raw phenotypes vs repeat length on '
         f"QC'ed sample subset, no covariates included, phenotype measured in {unit}"),
        ('pcausal', 'FINEMAP posterior probability causality'),
        ('previously_reported_association', 'Whether or not this row was included in the table '
         'because it was reported previously in the literature')
    ]

    readme.write(
        'table.tab contains one row for each variant that was not filtered prior to association '
        'testing and either was previously reported as an association or both had'
        'association p-value >= 5e-8 and'
        'FINEMAP posterior probability of causaility >= 0.05\n'
        '\n'
        'table.tab contains the following columns:\n'
    )
    for column, description in columns:
        readme.write(column)
        if description is not None:
            readme.write(' - ')
            readme.write(description)
        readme.write('\n')
    readme.write('\n')

    # chrom, signal_region, SNPSTR_start_pos, pcausal
    signals = pd.read_csv(
        f'{ukb}/finemapping/finemap_results/{phenotype}/summary/all_STR_contribs.tab',
        skiprows = 1,
        delimiter='\t'
    )
    rename_column_pd(signals, 'signal', 'signal_region')

    signal_split = signals['signal_region'].str.split('_', n=1, expand=True)
    signals['chrom'] = signal_split[0].astype(int)
    str_split = signals['STR'].str.split('_', n=1, expand=True)
    signals['pos'] = str_split[1].astype(int)

    str_results_fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'

    '''
    chrom
    SNPSTR_start_pos
    repeat_unit
    period
    alleles
    ref_len
    total_per_allele_dosages
    total_hardcall_alleles
    total_hardcall_genotypes
    subset_total_per_allele_dosages
    subset_total_hardcall_alleles
    subset_total_hardcall_genotypes
    subset_het
    subset_entropy
    subset_HWEP
    subset_allele_dosage_r2
    association_p_value
    '''
    associations = pd.read_csv(
        str_results_fname,
        header=0,
        delimiter='\t',
        dtype=utils.get_dtypes(str_results_fname, {'locus_filtered': str})
    )
    rename_column_pd(associations, 'ref_len', 'reference_allele')
    rename_column_pd(associations, 'motif', 'repeat_unit')
    rename_column_pd(associations, 'subset_het', 'subset_heterozygosity')
    rename_column_pd(associations, f'p_{phenotype}', 'association_p_value')
    for STR in previous_STRs:
        chrom, pos = STR.split(':')
        if not np.any((associations['chrom'] == int(chrom)) & (associations['pos'] == int(pos))):
            readme.write(f'Never ran association on previously reported STR association at {STR}\n')
            previous_STRs[STR] = True

    signals = signals.merge(
        associations,
        on=['chrom', 'pos'],
        how='right'
    )
    previous_STR_rows = None
    for STR, already_dropped in previous_STRs.items():
        if already_dropped:
            continue
        chrom, pos = STR.split(':')
        new_rows = (signals['chrom'] == int(chrom)) & (signals['pos'] == int(pos))
        if previous_STR_rows is None:
            previous_STR_rows = new_rows
        else:
            previous_STR_rows |= new_rows
    if previous_STR_rows is None:
        previous_STR_rows = np.full((signals.shape[0]), False, dtype=bool)

    signals['previously_reported_association'] = previous_STR_rows
    signals = signals[
        previous_STR_rows | (
            ~np.isnan(signals['pcausal'].to_numpy()) & (signals['pcausal'] >= 0.05)
    )]

    signals = signals[signals['association_p_value'] < 5e-8]
    signals.reset_index(inplace=True)
    nrows = signals.shape[0]

    for STR, already_dropped in previous_STRs.items():
        chrom, pos = STR.split(':')
        if not already_dropped and not np.any((signals['chrom'] == int(chrom)) & (signals['pos'] == int(pos))):
            readme.write(f'Omitting previously reported STR association at {STR} as it did not pass the GWAS association significance threshold\n')
            previous_STRs[STR] = True

    rename_column_pd(signals, 'pos', 'SNPSTR_start_pos')

    # calculate effect sizes and third largest allele freq
    dphen_units = np.full(nrows, np.nan)
    dphen_sds = np.full(nrows, np.nan)
    third_on_allele_freqs = np.full(nrows, np.nan)
    for idx in range(nrows):
        single_dosages = prep_dict(ast.literal_eval(signals['subset_total_per_allele_dosages'][idx]))
        single_dosages = prep_dict(ast.literal_eval(signals['subset_total_per_allele_dosages'][idx]))
        total_dosage = sum(single_dosages.values())
        single_freqs = {key: value/total_dosage for key, value in single_dosages.items()}

        if len(single_freqs) <= 2:
            third_on_allele_freqs[idx] = 0
        else:
            third_on_allele_freqs[idx] = sum(sorted(single_freqs.values())[:-2])

        doubled_freqs = {}
        for k1 in single_freqs:
            for k2 in single_freqs:
                key = k1 + k2
                if key not in doubled_freqs:
                    doubled_freqs[key] = 0
                doubled_freqs[key] += single_freqs[k1]*single_freqs[k2]

        phenotype_vals = prep_dict(ast.literal_eval(signals[f'mean_{phenotype}_per_single_dosage'][idx]))
        for key in phenotype_vals:
            if key not in doubled_freqs:
                print(phenotype_vals)
                print(doubled_freqs)
                assert False
        for key in doubled_freqs:
            if key not in phenotype_vals:
                if not doubled_freqs[key] <= 0.00001:
                    print(phenotype_vals)
                    print(doubled_freqs)
                    assert False

        xs = sorted(phenotype_vals.keys())
        ys = np.array([phenotype_vals[x] for x in xs])
        weights = np.array([doubled_freqs[x] for x in xs])
        xs = np.array(xs).reshape(-1, 1)
        xs = np.concatenate((xs, np.ones((xs.shape[0], 1))), axis=1)

        model = WLS(ys, xs, weights)
        reg_result = model.fit()
        coef = reg_result.params[0]
        scaled_coef = coef*DescrStatsW(xs[:, 0], weights=weights).std
        dphen_units[idx] = coef
        dphen_sds[idx] = scaled_coef
    signals['subset_multiallelicness'] = third_on_allele_freqs
    signals[col_dphen_unit] = dphen_units
    signals[col_dphen_sd] = dphen_sds

    chroms = []
    SNPSTR_poses = []
    SNPSTR_IDs = []
    starts = []
    ends = []
    for chrom in range(1, 23):
        vcf = cyvcf2.VCF(f'{ukb}/snpstr/info_field/chr{chrom}.vcf.gz')
        for variant in vcf:
            chroms.append(int(variant.CHROM))
            SNPSTR_poses.append(variant.POS)
            SNPSTR_IDs.append(variant.ID)
            starts.append(variant.INFO['START'])
            ends.append(variant.INFO['END'])
        vcf.close()

    info_df = pd.DataFrame.from_dict({
        'chrom': chroms,
        'SNPSTR_start_pos': SNPSTR_poses,
        'SNPSTR_ID': SNPSTR_IDs,
        'start_pos': starts,
        'end_pos': ends
    })

    signals = signals.merge(
        info_df,
        on=['chrom', 'SNPSTR_start_pos'],
        how="left"
    )

    signals = signals[[colname for colname, description in columns]]

    signals.to_csv(
        f'{ukb}/finemapping/summary/{phenotype}_table.tab',
        sep='\t',
        index=False,
        na_rep = 'NA'
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('--previous-STR-findings', nargs='+', default=[])
    args = parser.parse_args()
    phenotype = args.phenotype
    previous_STRs = {STR : False for STR in args.previous_STR_findings}

    with open(f'{ukb}/finemapping/summary/{phenotype}_table_README.txt', 'w') as readme:
        main(readme, phenotype, previous_STRs)
