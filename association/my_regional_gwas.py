#!/usr/bin/env python3

import argparse
import datetime
import os
import shutil
import tempfile
import time

import cyvcf2
import numpy as np
from statsmodels.regression.linear_model import OLS
import statsmodels.stats.weightstats

import load_and_filter_genotypes
import python_array_utils as utils

ukb = os.environ['UKB']

def perform_regional_gwas_helper(phenotype, outfile, runtype, imputation_run_name, region):
    outfile.write("chrom\tpos\talleles\tlocus_filtered\t"
                  f"p_{phenotype}\tcoeff_{phenotype}\tcoeff_intercept\tR^2\t")
    outfile.flush()

    n_loci = 0
    batch_time = 0
    batch_size = 50
    total_time = 0

    residuals = np.load(f'{ukb}/traits/adjusted_srin_phenotypes/{phenotype}_linear.npy')

    # order samples according to order in genetics files
    bgen_samples = []
    with open(f'{ukb}/microarray/ukb46122_hap_chr1_v2_s487314.sample') as samplefile:
        for num, line in enumerate(samplefile):
            if num <= 1:
                # skip first two lines
                continue
            bgen_samples.append(line.split()[0])
    assert len(bgen_samples) == 487409
    samples_array = np.array(bgen_samples, dtype=float).reshape(-1, 1)
    merge = utils.merge_arrays(samples_array, residuals)
    unfiltered_samples = ~np.isnan(merge[:, 1])

    residuals = merge[unfiltered_samples, 1]

    ori_phenotypes = np.load(f'{ukb}/traits/phenotypes/{phenotype}.npy')
    ori_phenotypes = utils.merge_arrays(samples_array, ori_phenotypes)[:, 1]
    ori_phenotypes = ori_phenotypes[unfiltered_samples]

    if runtype == 'strs':
        genotype_iter = load_and_filter_genotypes.load_strs(
            imputation_run_name, region, unfiltered_samples
        )
    elif runtype == 'imputed-snps':
        genotype_iter = load_and_filter_genotypes.load_imputed_snps(
            region, unfiltered_samples
        )
    else:
        raise ValueError("not implemented for this runtype")

    # first yield is special
    extra_detail_fields = next(genotype_iter)
    outfile.write('\t'.join(extra_detail_fields) + '\t')
    outfile.write(f'mean_residual_{phenotype}_per_single_dosage\t'
                  '0.05_significance_CI\t'
                  '5e-8_significance_CI')

    if runtype == 'strs':
        outfile.write(f'\tmean_residual_{phenotype}_per_paired_dosage\t'
                      '0.05_significance_CI\t'
                      '5e-8_significance_CI')
    outfile.write('\n')
    outfile.flush()

    start_time = time.time()
    for dosage_gts, unique_alleles, chrom, pos, locus_filtered, locus_details in genotype_iter:
        assert len(locus_details) == len(extra_detail_fields)

        n_loci += 1
        allele_names = ','.join(list(unique_alleles.astype(str)))
        outfile.write(f"{chrom}\t{pos}\t{allele_names}\t")
        if locus_filtered:
            outfile.write(f'{locus_filtered}\t1\tnan\tnan\tnan\t')
            outfile.write('\t'.join(locus_details))
            if runtype == 'strs':
                outfile.write('\tnan'*6 + '\n')
            else:
                outfile.write('\tnan'*3 + '\n')
            outfile.flush()
            continue
        else:
            outfile.write('False\t')

        if runtype == 'strs':
            gts = np.sum([_len*np.sum(dosages, axis=1) for
                          _len, dosages in dosage_gts.items()], axis=0)
        else:
            gts = dosage_gts[:, 1] + 2*dosage_gts[:, 2]
        gt_const = np.concatenate((
            gts.reshape(-1, 1),  np.ones((gts.shape[0], 1))
        ), axis = 1)

        #do da regression
        model = OLS(
            residuals,
            gt_const,
            missing='drop'
        )
        reg_result = model.fit()
        pval = reg_result.pvalues[0]
        coef = reg_result.params[0]
        intercept_coef = reg_result.params[1]
        rsquared = reg_result.rsquared
        outfile.write(f"{pval:.2e}\t{coef}\t{intercept_coef}\t{rsquared}\t")
        outfile.write('\t'.join(locus_details) + '\t')

        if runtype == 'strs':
            single_dosages = {}

            paired_dosages = {}
            for len1 in unique_alleles:
                for len2 in unique_alleles:
                    if len1 != len2:
                        dosages = (dosage_gts[len1][:, 0]*dosage_gts[len2][:, 1] +
                                   dosage_gts[len1][:, 1]*dosage_gts[len2][:, 0])
                    else:
                        dosages = dosage_gts[len1][:, 0]*dosage_gts[len1][:, 1]
                    if len1 + len2 not in single_dosages:
                        single_dosages[len1 + len2] = dosages
                    else:
                        single_dosages[len1 + len2] += dosages
                    minlen = min(len1, len2)
                    maxlen = max(len1, len2)
                    paired_dosages[(minlen, maxlen)] = dosages
            single_dosage_means = {}
            single_dosage_95_CI = {}
            single_dosage_GWAS_CI= {}
            for _len, dosages in single_dosages.items():
                if len(np.unique(ori_phenotypes[dosages != 0])) <= 1:
                    continue
                mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                    ori_phenotypes,
                    weights = dosages
                )
                single_dosage_means[_len] = mean_stats.mean
                single_dosage_95_CI[_len] = mean_stats.tconfint_mean()
                single_dosage_GWAS_CI[_len] = mean_stats.tconfint_mean(5e-8)
            paired_dosage_means = {}
            paired_dosage_95_CI = {}
            paired_dosage_GWAS_CI= {}
            for _len, dosages in paired_dosages.items():
                if len(np.unique(ori_phenotypes[dosages != 0])) <= 1:
                    continue
                mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                    ori_phenotypes,
                    weights = dosages
                )
                paired_dosage_means[_len] = mean_stats.mean
                paired_dosage_95_CI[_len] = mean_stats.tconfint_mean()
                paired_dosage_GWAS_CI[_len] = mean_stats.tconfint_mean(5e-8)
            outfile.write(load_and_filter_genotypes.dict_str(single_dosage_means) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(single_dosage_95_CI) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(single_dosage_GWAS_CI) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(paired_dosage_means) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(paired_dosage_95_CI) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(paired_dosage_GWAS_CI) + '\n')
        else:
            single_dosage_means = {}
            single_dosage_95_CI = {}
            single_dosage_GWAS_CI= {}
            for alt_count in range(3):
                mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                    ori_phenotypes,
                    weights = dosage_gts[:, alt_count]
                )
                single_dosage_means[alt_count] = mean_stats.mean
                single_dosage_95_CI[alt_count] = mean_stats.tconfint_mean()
                single_dosage_GWAS_CI[alt_count] = mean_stats.tconfint_mean(5e-8)
            outfile.write(load_and_filter_genotypes.dict_str(single_dosage_means) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(single_dosage_95_CI) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(single_dosage_GWAS_CI) + '\n')

        outfile.flush()

        duration = time.time() - start_time
        total_time += duration
        batch_time += duration
        if n_loci % batch_size == 0:
            print(
                f"time/locus (last {batch_size}): "
                f"{batch_time/batch_size}s\n"
                f"time/locus ({n_loci} total loci): {total_time/n_loci}s\n",
                flush = True
            )
            batch_time = 0
        start_time = time.time()
    if n_loci > 0:
        print(
            f"Done.\nTotal loci: {n_loci}\nTotal time: {total_time}s\ntime/locus: {total_time/n_loci}s\n",
            flush=True
        )
    else:
        print(f"No variants found in the region {region}\n", flush=True)


def perform_regional_gwas(phenotype, region, runtype, imputation_run_name):
    chrom, poses = region.split(':')
    start, end = poses.split('-')
    region_str = f'{chrom}_{start}_{end}'
    if runtype == 'strs':
        dirname = 'str'
    elif runtype == 'imputed-snps':
        dirname = 'imputed_snp'

    with tempfile.NamedTemporaryFile(mode='w+') as outfile:
        print(f"Writing output to temp file {outfile.name}", flush=True)
        perform_regional_gwas_helper(
            phenotype, outfile, runtype, imputation_run_name, region
        )
        shutil.copy(
            outfile.name,
            f'{ukb}/association/results/{phenotype}/my_{dirname}/batches/chr{region_str}.tab'
        )

def write_str_readme(phenotype, imputation_run_name):
    with open(f'{ukb}/association/results/{phenotype}/my_str/README.txt', 'w') as readme:
        today = datetime.datetime.now().strftime("%Y_%m_%d")
        readme.write(f"Run date: {today}\n")

        readme.write(f"Working with strs from imputation run {imputation_run_name}\n")
        readme.write("Working with length dosages, no call level filters.\n")
        readme.write("Filtering loci with fewer than 20 minor allele hardcalls.\n")
        readme.write("(Note: hardcalls mentioned in locus details are phased hardcalls,"
                     " and in rare cases will not correspond to the maximum likelihood "
                     "unphased allele)\n")
        readme.flush()

def write_imputed_snp_readme(phenotype):
    with open(f'{ukb}/association/results/{phenotype}/my_imputed_snp/README.txt', 'w') as readme:
        today = datetime.datetime.now().strftime("%Y_%m_%d")
        readme.write(f"Run date: {today}\n")

        readme.write("Working with dosages of alternate allele, no call level filters ")
        readme.write("(so allele 0 corresponds to reference, 1 to alternate).\n")
        readme.write("Filtering loci with fewer than 20 minor allele hardcalls.\n")
        readme.flush()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'runtype',
        choices=['imputed-snps', 'strs']
    )
    parser.add_argument('phenotype')
    parser.add_argument('--readme', action='store_true')
    parser.add_argument('--region')
    parser.add_argument('--imputation-run-name')
    args = parser.parse_args()

    assert args.readme == (args.region is None)
    assert (args.imputation_run_name is not None) == (args.runtype == 'strs')

    if args.readme:
        if args.runtype == 'strs':
            write_str_readme(args.phenotype, args.imputation_run_name)
        elif args.runtype == 'imputed-snps':
            write_imputed_snp_readme(args.phenotype)
        else:
            raise ValueError("Readme not implemented for this runtype")
        return
    else:
        perform_regional_gwas(
            args.phenotype, args.region, args.runtype, args.imputation_run_name
        )


if __name__ == '__main__':
    main()

