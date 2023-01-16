#!/usr/bin/env python3

import argparse
import datetime
import shutil
import tempfile
import time

import numpy as np
import polars as pl
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS
import statsmodels.stats.weightstats

import load_and_filter_genotypes
import python_array_utils as utils
import sample_utils
import weighted_binom_conf

def perform_regional_gwas_helper(
    outfile,
    pheno_and_covars_fname,
    shared_covars_fname,
    untransformed_phenotypes_fname,
    all_samples_fname,
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
    bgen_samples = sample_utils.get_samples(all_samples_fname)
    assert len(bgen_samples) == 487409
    samples_array = np.array(bgen_samples, dtype=float).reshape(-1, 1)
    merge = utils.merge_arrays(samples_array, covars)
    unfiltered_samples = ~np.any(np.isnan(merge), axis=1)

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

    # TODO maybe only do these calculations if there's a flag on
    if runtype == 'strs':
        outfile.write(
            'total_subset_dosage_per_summed_gt\t'
            f'{stat}_{phenotype}_per_summed_gt\t'
            'summed_0.05_significance_CI\t'
            'summed_5e-8_significance_CI\t'
            f'{stat}_{phenotype}_residual_per_summed_gt\t'
            'res_per_sum_0.05_significance_CI\t'
            'res_per_sum_5e-8_significance_CI\t'
            'total_subset_dosage_per_paired_gt\t'
            f'{stat}_{phenotype}_per_paired_gt\t'
            'paired_0.05_significance_CI\t'
            'paired_5e-8_significance_CI\t'
            f'{stat}_{phenotype}_residual_per_paired_gt\t'
            'res_per_paired_0.05_significance_CI\t'
            'res_per_paired_5e-8_significance_CI'
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
            gts = np.sum([len_*np.sum(dosages, axis=1) for
                          len_, dosages in dosage_gts.items()], axis=0)
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
            outfile.write(f"{pval:.2e}\t{coef/std}\t{se/std}\t{rsquared}\t")
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

        outfile.write('\t'.join(locus_details))

        if runtype == 'strs':
            dosages_per_summed_gt = {}

            dosages_per_paired_gt = {}
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
                    summedlen_ = round(len1 + len2, 2)
                    if summedlen_ not in dosages_per_summed_gt:
                        dosages_per_summed_gt[summedlen_] = dosages
                    else:
                        dosages_per_summed_gt[summedlen_] += dosages
                    minlen = min(len1, len2)
                    maxlen = max(len1, len2)
                    dosages_per_paired_gt[(minlen, maxlen)] = dosages

            outfile.write('\t' + load_and_filter_genotypes.dict_str({key: np.sum(arr) for key, arr in dosages_per_summed_gt.items()}))
            if not binary or binary == 'linear':
                #do da regression
                untrans_model = OLS(
                    ori_phenotypes,
                    covars[:, 1:],
                    missing='drop',
                )
            else:
                untrans_model = sm.GLM(
                    ori_phenotypes,
                    covars[:, 1:],
                    missing='drop',
                    family=sm.families.Binomial()
                )
            untrans_reg_result = untrans_model.fit()
            residual_phenotypes = ori_phenotypes - untrans_reg_result.fittedvalues

            for phenotypes in ori_phenotypes, residual_phenotypes:
                summed_gt_stat = {}
                summed_gt_95_CI = {}
                summed_gt_GWAS_CI= {}
                if not binary:
                    for len_, dosages in dosages_per_summed_gt.items():
                        if len(np.unique(phenotypes[dosages != 0])) <= 1:
                            continue
                        mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                            phenotypes,
                            weights = dosages
                        )
                        summed_gt_stat[len_] = mean_stats.mean
                        summed_gt_95_CI[len_] = mean_stats.tconfint_mean()
                        summed_gt_GWAS_CI[len_] = mean_stats.tconfint_mean(5e-8)
                else:
                    for len_, dosages in dosages_per_summed_gt.items():
                        if not np.any(dosages != 0):
                            continue
                        p, lower, upper = weighted_binom_conf.weighted_binom_conf(
                            dosages, phenotypes, 0.05
                        )
                        summed_gt_stat[len_] = p
                        summed_gt_95_CI[len_] = (lower, upper)
                        _,  lower_gwas, upper_gwas = weighted_binom_conf.weighted_binom_conf(
                            dosages, phenotypes, 5e-8
                        )
                        summed_gt_GWAS_CI[len_] = (lower_gwas, upper_gwas)
                outfile.write('\t' + load_and_filter_genotypes.dict_str(summed_gt_stat))
                outfile.write('\t' + load_and_filter_genotypes.dict_str(summed_gt_95_CI))
                outfile.write('\t' + load_and_filter_genotypes.dict_str(summed_gt_GWAS_CI))

            outfile.write('\t' + load_and_filter_genotypes.dict_str({key: np.sum(arr) for key, arr in dosages_per_paired_gt.items()}))
            for phenotypes in ori_phenotypes, residual_phenotypes:
                paired_gt_stat = {}
                paired_gt_95_CI = {}
                paired_gt_GWAS_CI= {}
                if not binary:
                    for len_, dosages in dosages_per_paired_gt.items():
                        if len(np.unique(phenotypes[dosages != 0])) <= 1:
                            continue
                        mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                            phenotypes,
                            weights = dosages
                        )
                        paired_gt_stat[len_] = mean_stats.mean
                        paired_gt_95_CI[len_] = mean_stats.tconfint_mean()
                        paired_gt_GWAS_CI[len_] = mean_stats.tconfint_mean(5e-8)
                else:
                    for len_, dosages in dosages_per_paired_gt.items():
                        if not np.any(dosages != 0):
                            continue
                        p, lower, upper = weighted_binom_conf.weighted_binom_conf(
                            dosages, phenotypes, 0.05
                        )
                        paired_gt_stat[len_] = p
                        paired_gt_95_CI[len_] = (lower, upper)
                        _,  lower_gwas, upper_gwas = weighted_binom_conf.weighted_binom_conf(
                            dosages, phenotypes, 5e-8
                        )
                        paired_gt_GWAS_CI[len_] = (lower_gwas, upper_gwas)
                outfile.write('\t' + load_and_filter_genotypes.dict_str(paired_gt_stat))
                outfile.write('\t' + load_and_filter_genotypes.dict_str(paired_gt_95_CI))
                outfile.write('\t' + load_and_filter_genotypes.dict_str(paired_gt_GWAS_CI))

        outfile.write('\n')
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

def get_genotype_iter_vars_file(str_vcf, vars_fname, samples):
    with open(vars_fname) as vars_file:
        next(vars_file)
        try:
            next(vars_file)
        except StopIteration:
            # this is a empty vars file
            # yield only the list of details fields and then exit without yielding variants
            itr = load_and_filter_genotypes.load_strs(
                str_vcf, f'1:1-1', samples
            )
            yield next(itr)
            return
    f = pl.read_csv(vars_file, sep='\t')
    chroms = f['chrom']
    poses = f['pos']
    first = True
    for (chrom, pos) in zip(chroms, poses):
        itr = load_and_filter_genotypes.load_strs(
            str_vcf, f'{chrom}:{pos}-{pos}', samples
        )
        # yield or skip the extra details line
        if first:
            yield next(itr)
            first = False
        else:
            next(itr)
        # yield the genotype
        yield next(itr)


def perform_regional_gwas(outfile, pheno_and_covars_fname, shared_covars_fname, untransformed_phenotypes_fname, all_samples_fname, phenotype, binary, region, vars_file, runtype, str_vcf, snp_bgen, snp_mfi, conditional_covars_fname, temp_dir):
    if runtype == 'strs':
        if region is not None:
            get_genotype_iter = lambda samples: load_and_filter_genotypes.load_strs(
                str_vcf, region, samples
            )
        else:
            assert vars_file is not None
            get_genotype_iter = lambda samples: get_genotype_iter_vars_file(str_vcf, vars_file, samples)
    elif runtype == 'imputed-snps':
        assert vars_file is None
        get_genotype_iter = lambda samples: load_and_filter_genotypes.load_imputed_snps(
            snp_bgen, snp_mfi, region, samples
        )
    else:
        raise ValueError("not implemented for this runtype")

    with tempfile.NamedTemporaryFile(dir=temp_dir, mode='w+') as temp_outfile:
        print(f"Writing output to temp file {temp_outfile.name}", flush=True)
        perform_regional_gwas_helper(
            temp_outfile, pheno_and_covars_fname, shared_covars_fname, untransformed_phenotypes_fname, all_samples_fname, get_genotype_iter, phenotype, binary, region, runtype, conditional_covars_fname
        )

        print(f"Copying {temp_outfile.name} to {outfile}")
        shutil.copy(
            temp_outfile.name,
            outfile
        )
        print("Done.")

def write_str_readme(outfile, str_vcf, binary):
    with open(outfile, 'w') as readme:
        today = datetime.datetime.now().strftime("%Y_%m_%d")
        readme.write(f"Run date: {today}\n")

        readme.write(f"Working with strs from imputation run {str_vcf}\n")
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
    parser.add_argument('--readme', action='store_true', default=False)
    parser.add_argument('--temp-dir')
    parser.add_argument('--region')
    parser.add_argument('--vars-file')
    parser.add_argument('--pheno-and-covars')
    parser.add_argument('--shared-covars')
    parser.add_argument('--conditional-covars')
    parser.add_argument('--untransformed-phenotypes')
    parser.add_argument('--all-samples-fname')
    parser.add_argument('--str-vcf')
    parser.add_argument('--snp-bgen')
    parser.add_argument('--snp-mfi')
    parser.add_argument('--binary', default=False, choices={'linear', 'logistic'})
    args = parser.parse_args()

    assert args.readme == (args.region is None and args.vars_file is None)
    assert (args.str_vcf is not None) == (args.runtype == 'strs')
    assert (not args.readme) or ((args.snp_bgen is not None) == (args.snp_mfi is not None) == (args.runtype == 'imputed-snps'))

    assert (args.region is not None) + (args.vars_file is not None) <= 1

    # this would require firth regression which I've removed
    assert not ((args.runtype == 'imputed-snps') and (args.binary == 'logistic'))

    assert (
        args.readme ==
        (args.pheno_and_covars is None) ==
        (args.shared_covars is None) ==
        (args.untransformed_phenotypes is None) ==
        (args.all_samples_fname is None) ==
        (args.temp_dir is None)
    )

    if args.conditional_covars is not None:
        assert args.region is not None

    if args.readme:
        if args.runtype == 'strs':
            write_str_readme(args.outfile, args.str_vcf, args.binary)
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
            args.all_samples_fname,
            args.phenotype,
            args.binary,
            args.region,
            args.vars_file,
            args.runtype,
            args.str_vcf,
            args.snp_bgen,
            args.snp_mfi,
            args.conditional_covars,
            args.temp_dir
        )


if __name__ == '__main__':
    main()

