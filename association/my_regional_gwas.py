#!/usr/bin/env python3

import argparse
import datetime
import os
import shutil
import tempfile
import time

import numpy as np
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS
import statsmodels.stats.weightstats

import load_and_filter_genotypes
import python_array_utils as utils
import sample_utils
import weighted_binom_conf

project_temp = os.environ['PROJECT_TEMP']

def perform_regional_gwas_helper(
    outfile,
    pheno_and_covars_fname,
    shared_covars_fname,
    untransformed_phenotypes_fname,
    get_genotype_iter,
    phenotype,
    binary,
    region,
    runtype,
    conditional_covars_fname = None
):

    outfile.write("chrom\tpos\talleles\tlocus_filtered\t"
                  f"p_{phenotype}\tcoeff_{phenotype}\t")
    if binary != 'logistic':
        outfile.write(f'se_{phenotype}\tR^2\t')
    else:
        outfile.write("unused_col\tunused_col\t")
    outfile.flush()

    n_loci = 0
    batch_time = 0
    batch_size = 50
    total_time = 0

    pheno_specific_covars = np.load(pheno_and_covars_fname)
    shared_covars = np.load(shared_covars_fname)
    covars = utils.merge_arrays(pheno_specific_covars, shared_covars)

    if conditional_covars_fname:
        gt_covars = np.load(conditional_covars_fname)
        covars = utils.merge_arrays(covars, gt_covars)

    # order samples according to order in genetics files
    bgen_samples = sample_utils.get_all_samples()
    assert len(bgen_samples) == 487409
    samples_array = np.array(bgen_samples, dtype=float).reshape(-1, 1)
    merge = utils.merge_arrays(samples_array, covars)
    unfiltered_samples = ~np.isnan(merge[:, 1])

    outcome = merge[unfiltered_samples, 1].copy()
    covars = merge[unfiltered_samples, :]
    covars = (covars - np.mean(covars, axis=0))/np.std(covars, axis=0)
    covars[:, 1] = 1 # reuse the column that was the outcome as the intercept
    
    ori_phenotypes = np.load(untransformed_phenotypes_fname)
    ori_phenotypes = utils.merge_arrays(samples_array, ori_phenotypes)[:, 1]
    ori_phenotypes = ori_phenotypes[unfiltered_samples]

    # first yield is special
    genotype_iter = get_genotype_iter(unfiltered_samples)
    extra_detail_fields = next(genotype_iter)
    outfile.write('\t'.join(extra_detail_fields) + '\t')

    if not binary:
        stat = 'mean'
    else:
        stat = 'fraction'

    outfile.write(f'{stat}_{phenotype}_per_single_dosage\t'
                  '0.05_significance_CI\t'
                  '5e-8_significance_CI')

    if runtype == 'strs':
        outfile.write(
            '\ttotal_subset_dosage_per_summed_gt\t'
            f'{stat}_{phenotype}_per_paired_dosage\t'
            '0.05_significance_CI\t'
            '5e-8_significance_CI'
        )
    outfile.write('\n')
    outfile.flush()

    start_time = time.time()
    for dosage_gts, unique_alleles, chrom, pos, locus_filtered, locus_details in genotype_iter:
        assert len(locus_details) == len(extra_detail_fields)

        covars[:, 0] = np.nan # reuse the column that was the ids as the genotypes

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
        std = np.std(gts)
        gts = (gts - np.mean(gts))/np.std(gts)
        covars[:, 0] = gts

        if not binary or binary == 'linear':
            #do da regression
            model = OLS(
                outcome,
                covars,
                missing='drop',
            )
            reg_result = model.fit()
            pval = reg_result.pvalues[0]
            coef = reg_result.params[0]
            se = reg_result.bse[0]
            rsquared = reg_result.rsquared
            outfile.write(f"{pval:.2e}\t{coef/std}\t{se}\t{rsquared}\t")
        else:
            model = sm.GLM(
                outcome,
                covars,
                missing='drop',
                family=sm.families.Binomial()
            )
            reg_result = model.fit()
            pval = reg_result.pvalues[0]
            coef = reg_result.params[0]
            outfile.write(f'{pval:.2e}\t{coef/std}\tnan\tnan\t')

        outfile.write('\t'.join(locus_details) + '\t')

        if runtype == 'strs':
            single_dosages = {}

            paired_dosages = {}
            for len1 in unique_alleles:
                for len2 in unique_alleles:
                    if len1 > len2:
                        continue
                    if len1 != len2:
                        dosages = (dosage_gts[len1][:, 0]*dosage_gts[len2][:, 1] +
                                   dosage_gts[len1][:, 1]*dosage_gts[len2][:, 0])
                    else:
                        dosages = dosage_gts[len1][:, 0]*dosage_gts[len1][:, 1]
                    if np.sum(dosages) <= 0:
                        continue
                    if len1 + len2 not in single_dosages:
                        single_dosages[len1 + len2] = dosages
                    else:
                        single_dosages[len1 + len2] += dosages
                    minlen = min(len1, len2)
                    maxlen = max(len1, len2)
                    paired_dosages[(minlen, maxlen)] = dosages
            single_dosage_stat = {}
            single_dosage_95_CI = {}
            single_dosage_GWAS_CI= {}
            paired_dosage_stat = {}
            paired_dosage_95_CI = {}
            paired_dosage_GWAS_CI= {}
            if not binary:
                for _len, dosages in single_dosages.items():
                    if len(np.unique(ori_phenotypes[dosages != 0])) <= 1:
                        continue
                    mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                        ori_phenotypes,
                        weights = dosages
                    )
                    single_dosage_stat[_len] = mean_stats.mean
                    single_dosage_95_CI[_len] = mean_stats.tconfint_mean()
                    single_dosage_GWAS_CI[_len] = mean_stats.tconfint_mean(5e-8)
                for _len, dosages in paired_dosages.items():
                    if len(np.unique(ori_phenotypes[dosages != 0])) <= 1:
                        continue
                    mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                        ori_phenotypes,
                        weights = dosages
                    )
                    paired_dosage_stat[_len] = mean_stats.mean
                    paired_dosage_95_CI[_len] = mean_stats.tconfint_mean()
                    paired_dosage_GWAS_CI[_len] = mean_stats.tconfint_mean(5e-8)
            else:
                for _len, dosages in single_dosages.items():
                    if not np.any(dosages != 0):
                        continue
                    p, lower, upper = weighted_binom_conf.weighted_binom_conf(
                        dosages, ori_phenotypes, 0.05
                    )
                    single_dosage_stat[_len] = p
                    single_dosage_95_CI[_len] = (lower, upper)
                    _,  lower_gwas, upper_gwas = weighted_binom_conf.weighted_binom_conf(
                        dosages, ori_phenotypes, 5e-8
                    )
                    single_dosage_GWAS_CI[_len] = (lower_gwas, upper_gwas)
                for _len, dosages in paired_dosages.items():
                    if not np.any(dosages != 0):
                        continue
                    p, lower, upper = weighted_binom_conf.weighted_binom_conf(
                        dosages, ori_phenotypes, 0.05
                    )
                    paired_dosage_stat[_len] = p
                    paired_dosage_95_CI[_len] = (lower, upper)
                    _,  lower_gwas, upper_gwas = weighted_binom_conf.weighted_binom_conf(
                        dosages, ori_phenotypes, 5e-8
                    )
                    paired_dosage_GWAS_CI[_len] = (lower_gwas, upper_gwas)
            outfile.write(load_and_filter_genotypes.dict_str(single_dosage_stat) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(single_dosage_95_CI) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(single_dosage_GWAS_CI) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str({key: np.sum(arr) for key, arr in single_dosages.items()}) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(paired_dosage_stat) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(paired_dosage_95_CI) + '\t')
            outfile.write(load_and_filter_genotypes.dict_str(paired_dosage_GWAS_CI) + '\n')
        else:
            single_dosage_stat = {}
            single_dosage_95_CI = {}
            single_dosage_GWAS_CI= {}
            if not binary:
                for alt_count in range(3):
                    mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                        ori_phenotypes,
                        weights = dosage_gts[:, alt_count]
                    )
                    single_dosage_stat[alt_count] = mean_stats.mean
                    single_dosage_95_CI[alt_count] = mean_stats.tconfint_mean()
                    single_dosage_GWAS_CI[alt_count] = mean_stats.tconfint_mean(5e-8)
            else:
                for alt_count in range(3):
                    p, lower, upper = weighted_binom_conf.weighted_binom_conf(
                        dosage_gts[:, alt_count], ori_phenotypes, 0.05
                    )
                    single_dosage_stat[alt_count] = p
                    single_dosage_95_CI[alt_count] = (lower, upper)
                    _,  lower_gwas, upper_gwas = weighted_binom_conf.weighted_binom_conf(
                        dosage_gts[:, alt_count], ori_phenotypes, 5e-8
                    )
                    single_dosage_GWAS_CI[alt_count] = (lower_gwas, upper_gwas)
            outfile.write(load_and_filter_genotypes.dict_str(single_dosage_stat) + '\t')
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


def perform_regional_gwas(outfile, pheno_and_covars_fname, shared_covars_fname, untransformed_phenotypes_fname, phenotype, binary, region, runtype, imputation_run_name, conditional_covars_fname):
    if runtype == 'strs':
        get_genotype_iter = lambda samples: load_and_filter_genotypes.load_strs(
            imputation_run_name, region, samples
        )
    elif runtype == 'imputed-snps':
        get_genotype_iter = lambda samples: load_and_filter_genotypes.load_imputed_snps(
            region, samples
        )
    else:
        raise ValueError("not implemented for this runtype")

    with tempfile.NamedTemporaryFile(dir=project_temp, mode='w+') as temp_outfile:
        print(f"Writing output to temp file {temp_outfile.name}", flush=True)
        perform_regional_gwas_helper(
            temp_outfile, pheno_and_covars_fname, shared_covars_fname, untransformed_phenotypes_fname, get_genotype_iter, phenotype, binary, region, runtype, conditional_covars_fname
        )

        print(f"Copying {temp_outfile.name} to {outfile}")
        shutil.copy(
            temp_outfile.name,
            outfile
        )
        print("Done.")

def write_str_readme(outfile, imputation_run_name, binary):
    with open(outfile, 'w') as readme:
        today = datetime.datetime.now().strftime("%Y_%m_%d")
        readme.write(f"Run date: {today}\n")

        readme.write(f"Working with strs from imputation run {imputation_run_name}\n")
        readme.write("Working with length dosages, no call level filters.\n")
        readme.write("Filtering loci with total non-major allele doesage less than 20.\n")

        if not binary:
            readme.write('Doing linear regressions against the continuous phenotype\n')
        elif binary == 'linear':
            readme.write('Doing linear regressions against the binary phenotype\n')
        else:
            readme.write('Doing logistic regressions against the binary phenotype. No longer '
                         'using firth penalized logistic regression when MAC <= 400, should but '
                         "this doesn't apply to any strs in this dataset. Instead, always using "
                         'standard logistic regression.\n')

def write_imputed_snp_readme(outfile, binary):
    with open(outfile, 'w') as readme:
        today = datetime.datetime.now().strftime("%Y_%m_%d")
        readme.write(f"Run date: {today}\n")

        readme.write("Working with dosages of alternate allele, no call level filters ")
        readme.write("(so allele 0 corresponds to reference, 1 to alternate).\n")
        readme.write("Filtering loci with total minor allele dosage less than 20.\n")

        if not binary:
            readme.write('Doing linear regressions against the continuous phenotype\n')
        elif binary == 'linear':
            readme.write('Doing linear regressions against the binary phenotype\n')
        else:
            readme.write('Doing logistic regressions against the binary phenotype. Using '
                         'Frith penalized logistic regression when MAC <= 400, '
                         'standard logistic regression otherwise.\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile')
    parser.add_argument(
        'runtype',
        choices=['imputed-snps', 'strs']
    )
    parser.add_argument('phenotype')
    parser.add_argument('--readme', action='store_true')
    parser.add_argument('--region')
    parser.add_argument('--pheno-and-covars')
    parser.add_argument('--shared-covars')
    parser.add_argument('--conditional-covars')
    parser.add_argument('--untransformed-phenotypes')
    parser.add_argument('--imputation-run-name')
    parser.add_argument('--binary', default=False, choices={'linear', 'logistic'})
    args = parser.parse_args()

    assert args.readme == (args.region is None)
    assert (args.imputation_run_name is not None) == (args.runtype == 'strs')

    # this would require firth regression which I've removed
    assert not ((args.runtype == 'imputed-snps') and (args.binary == 'logistic'))

    assert (
        (args.region is not None) ==
        (args.pheno_and_covars is not None) ==
        (args.shared_covars is not None) ==
        (args.untransformed_phenotypes is not None)
    )

    if args.conditional_covars is not None:
        assert args.region is not None

    if args.readme:
        if args.runtype == 'strs':
            write_str_readme(args.outfile, args.imputation_run_name, args.binary)
        elif args.runtype == 'imputed-snps':
            write_imputed_snp_readme(args.outfile, args.binary)
        else:
            raise ValueError("Readme not implemented for this runtype")
        return
    else:
        perform_regional_gwas(
            args.outfile,
            args.pheno_and_covars,
            args.shared_covars,
            args.untransformed_phenotypes,
            args.phenotype,
            args.binary,
            args.region,
            args.runtype,
            args.imputation_run_name,
            args.conditional_covars
        )


if __name__ == '__main__':
    main()

