#!/usr/bin/env python3

import argparse
import os
import sys

import cyvcf2
import numpy as np
import polars as pl
from statsmodels.regression.linear_model import OLS

import python_array_utils
import sample_utils

ukb = os.environ['UKB']

sys.path.insert(0, f'{ukb}/../trtools/repo')

import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

def var_dosage_gts(var, samp_idx, alleles):
    gts = np.zeros(np.sum(samp_idx))
    for allele in alleles:
        if allele > 0:
            gts += (var.format('AP1') + var.format('AP2'))[samp_idx, allele-1]
        else:
            assert allele == 0
            gts += (
                np.maximum(0, 1 - np.sum(var.format('AP1'), axis=1)) +
                np.maximum(0, 1 - np.sum(var.format('AP2'), axis=1))
            )[samp_idx]
    return gts

def str_dosage_gts(var, samp_idx):
    lens = [len(var.REF)] + [len(alt) for alt in var.ALT]
    gts = np.zeros(np.sum(samp_idx))
    gts += lens[0]*(
        np.maximum(0, 1 - np.sum(var.format('AP1'), axis=1)) +
        np.maximum(0, 1 - np.sum(var.format('AP2'), axis=1))
    )[samp_idx]
    for i in range(1, len(lens)):
        gts += lens[i]*(var.format('AP1') + var.format('AP2'))[samp_idx, i-1]
    return gts

def standardize(gts):
    if np.std(gts) == 0:
        return gts
    return (gts - np.mean(gts))/np.std(gts)

def main():
    # do an association test of the combined dosage of the listed alleles vs the listed phenotypes in each ethnicity
    parser = argparse.ArgumentParser()
    #parser.add_argument('chrom', type=int)
    #parser.add_argument('pos', type=int)
    parser.add_argument('var_file')
    #parser.add_argument('--phenotypes', nargs='+')
    #parser.add_argument('--alleles', type=int, nargs='+') #allele indicies, i.e. ths SNP is present in alleles 0 (ref), 2 (2nd alt) and 5
    args = parser.parse_args()

    print('str\tvariant\tethnicity\tvar frequency\tstr var r2\tphenotype\tstr p-val\tvar p-val\tstr p-val conditioning on var')

    scovs = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')

    variants = pl.read_csv(args.var_file, sep='\t')
    for i in range(variants.shape[0]):
        chrom = variants[i, 'chrom']
        pos = variants[i, 'pos']
        str_ = f'{chrom}:{pos}'
        var = next(cyvcf2.VCF(f'{ukb}/str_imputed/runs/first_pass/vcfs/annotated_strs/chr{chrom}.vcf.gz')(str_))
        alleles = [int(num) for num in variants[i, 'alleles'].split(',')]
        name = variants[i, 'name']
        phenotypes = variants[i, 'phenos'].split(',')
        for phenotype in phenotypes:
            for ethnicity in ('white_brits', 'black', 'south_asian', 'chinese', 'irish', 'white_other'):
                #ethnicity
                total_samp_idx = sample_utils.get_samples_idx_ethnicity(ethnicity)
                total_var_gts = var_dosage_gts(var, total_samp_idx, alleles)
                total_str_gts = str_dosage_gts(var, total_samp_idx)
                var_freq = np.sum(total_var_gts)/(2*total_var_gts.shape[0])
                corr = np.corrcoef(total_var_gts, total_str_gts)
                assert corr.shape == (2,2)
                str_var_r2 = corr[0,1]**2
                #ethnicity,pheno
                pcovs = np.load(f'{ukb}/traits/subset_transformed_phenotypes/{ethnicity}/{phenotype}.npy')
                samps = sample_utils.get_ordered_samples_phenotype(ethnicity, phenotype).reshape(-1, 1)
                covs = python_array_utils.merge_arrays(python_array_utils.merge_arrays(samps, pcovs), scovs)

                outcomes = covs[:, 1]
                covs = covs[:, 2:]

                samp_idx = sample_utils.get_samples_idx_phenotype(ethnicity, phenotype)
                var_gts = standardize(var_dosage_gts(var, samp_idx, alleles))
                str_gts = standardize(str_dosage_gts(var, samp_idx))
                str_best_guess_gts = trh.HarmonizeRecord(
                    vcfrecord=var, vcftype='beagle-hipstr'
                ).GetGenotypeIndicies()[samp_idx, :-1]

                str_p = OLS(
                    outcomes,
                    np.hstack((covs, np.ones((covs.shape[0], 1)), str_gts.reshape(-1, 1)))
                ).fit().pvalues[-1]

                if np.all(var_gts == 0) or np.all(var_gts == 2):
                    var_p = 1
                    str_cond_p = str_p
                else:
                    var_p = OLS(
                        outcomes,
                        np.hstack((covs, np.ones((covs.shape[0], 1)), var_gts.reshape(-1, 1)))
                    ).fit().pvalues[-1]
                    str_cond_p = OLS(
                        outcomes,
                        np.hstack((covs, np.ones((covs.shape[0], 1)), var_gts.reshape(-1, 1), str_gts.reshape(-1, 1)))
                    ).fit().pvalues[-1]

                print(f'{str_}\t{name}\t{ethnicity}\t{var_freq:.3g}\t{str_var_r2:.3g}\t{phenotype}\t{str_p:.3g}\t{var_p:.3g}\t{str_cond_p:.3g}')

if __name__ == '__main__':
    main()
# specify SNPs/indels by specific alleles
# frequency of those varaints, per population
# correlation between length and presence of those variants, per population
# p-values of association of those variants, per population and phenotype
# p-values of association of length with phenotype conditioned on the presence of those variants one at a time, per population and phenotype
# Be able to do the above where instead of conditioning on the presence of a variant, run two separate correlations - among those with covariate and those without
# be able to stratify plot of length vs phenotype association against presence of variant or no
# finemapping with variants included?
