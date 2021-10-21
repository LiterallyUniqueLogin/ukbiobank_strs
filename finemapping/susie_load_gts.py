#!/usr/bin/env python3

import argparse
import datetime
import math
import os
import pathlib
import subprocess as sp
import sys
import time

import numpy as np
import sklearn.linear_model

import load_and_filter_genotypes as lfg
import python_array_utils as utils

ukb = os.environ['UKB']

def get_str_dosages(str_dosages_dict):
    return np.sum([_len*np.sum(dosages, axis=1) for
                   _len, dosages in str_dosages_dict.items()], axis=0)

def load_gts(imputation_run_name, phenotype, chrom, start_pos, end_pos):
    filter_set_fname = f'{ukb}/finemapping/str_imp_snp_overlaps/chr{chrom}_to_filter.tab'
    outdir = f'{ukb}/finemapping/susie_results/{phenotype}/{chrom}_{start_pos}_{end_pos}'

    p_cutoff = 5e-4
   
    today = datetime.datetime.now().strftime("%Y_%M_%D")
    with open(f"{outdir}/README.txt", 'w') as readme:
        readme.write(
            f'Run date: {today}\n'
            'Loading STR and SNP gts and the phenotype values, regressing out '
            'covariates from all of those, then running SuSiE. '
            f'Variants for which association tests were skipped or with p >= {p_cutoff} are excluded. '
            'SNPs in the filter set are also skipped. '
            f'(Filter set at {filter_set_fname})\n'
        )

    with open(f'{ukb}/sample_qc/runs/{phenotype}/combined_unrelated.sample') as samples_file:
        samples = np.array([line.strip() for line in samples_file][1:], dtype=int).reshape(-1, 1)

    imp_snp_samples_filepath = f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample'
    with open(imp_snp_samples_filepath) as imp_snp_samples_file:
        imp_snp_samples = np.array([line.split()[0] for line in imp_snp_samples_file][2:], dtype=int).reshape(-1, 1)

    samples_indicator = np.concatenate((samples, samples), axis=1)
    samples_merge = utils.merge_arrays(imp_snp_samples, samples_indicator)
    assert samples_merge.shape[1] == 2
    sample_idx = ~np.isnan(samples_merge[:, 1])

    pheno_covars = np.load(f'{ukb}/traits/subset_transformed_phenotypes/{phenotype}.npy')
    shared_covars = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')
    covars = utils.merge_arrays(pheno_covars, shared_covars)

    # reorder samples to the proper order
    covars = utils.merge_arrays(samples_merge[sample_idx, 0:1], covars)
    pheno_vals = covars[:, 1]
    covars = covars[:, 2:]

    # first choose STRs and SNPs only with p < p_cutoff to lessen memory burden
    strs_to_include = set()
    snps_to_include = set()

    str_results_fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    with open(str_results_fname) as str_results:
        header = next(str_results).strip().split('\t')
    str_p_idx = header.index(f'p_{phenotype}')

    out = sp.run(
        f"grep -P '^{chrom}\t' {str_results_fname} | "
        f"awk '{{ if ({start_pos} <= $2 && $2 <= {end_pos} && ${str_p_idx + 1} < {p_cutoff}) "
        " { print $2 } }' ",
        shell=True,
        check = True,
        capture_output=True
    )
    for line in out.stdout.decode().split('\n'):
        strs_to_include.add(line.strip())
    if len(strs_to_include) == 0:
        with open(f"{outdir}/README.txt", 'w') as readme:
            readme.write(
                f"No STRs were found in the region with p < {p_cutoff}, "
                "so finemapping is being skipped.\n"
            )
        pathlib.Path(f'{outdir}/no_strs').touch()
        pathlib.Path(f'{outdir}/pheno_residuals.npy').touch()
        pathlib.Path(f'{outdir}/gt_residuals.npy').touch()
        sys.exit()

    snps_to_filter = set()
    with open(filter_set_fname) as filter_file:
        next(filter_file) # skip header
        for line in filter_file:
            pos, ref, alt = line.strip().split('\t')[3:6]
            snps_to_filter.add((pos, ref, alt))

    snp_results_fname = f'{ukb}/association/results/{phenotype}/plink_snp/results.tab'
    with open(snp_results_fname) as snp_results:
        header = next(snp_results).strip().split('\t')
    snp_p_idx = header.index('P')

    out = sp.run(
        f"grep -P '^{chrom}\t' {snp_results_fname} | "
        f"awk '{{ if ({start_pos} <= $2 && $2 <= {end_pos} && ${snp_p_idx + 1} != \"NA\" && ${snp_p_idx + 1} < {p_cutoff}) "
        " { print $2 \" \" $4 \" \" $5 } }' ",
        shell=True,
        check = True,
        capture_output=True
    )
    for line in out.stdout.decode().split('\n'):
        snp = tuple(line.strip().split())
        if snp not in snps_to_filter:
            snps_to_include.add(snp)

    var_names = []
    gts = []

    region = f'{chrom}:{start_pos}-{end_pos}'
    str_itr = lfg.load_strs(imputation_run_name, region, sample_idx)
    next(str_itr) # skip the locus details list

    print('loading STRs... ', flush=True)
    start = time.time()
    for str_count, (str_dosages_dict, _, _, str_pos, str_locus_filtered, _) in enumerate(str_itr):
        str_count += 1
        print(f'loading STR {str_count} ... ', flush=True)
        if str_locus_filtered:
            continue
        if str(str_pos) not in strs_to_include:
            continue
        var_names.append(f'STR_{str_pos}')
        gts.append(get_str_dosages(str_dosages_dict))
    print(f"Time: {time.time() - start}s")
    if len(gts[0]) != len(strs_to_include):
        print(var_names)
        print(strs_to_include)
        assert False

    snp_itr = lfg.load_imputed_snps(region, sample_idx)

    print('loading SNPs... ', flush=True)
    start = time.time()
    for snp_count, (snp_dosages, alleles, _, snp_pos, snp_filtered, _) in enumerate(snp_itr):
        snp_count += 1
        if snp_count % 10 == 0:
            print(f'loading SNP {snp_count} ... ', flush=True)
        if snp_filtered:
            continue
        if (str(snp_pos), *alleles) not in snps_to_include:
            continue
        var_names.append(f'SNP_{snp_pos}_' + '_'.join(alleles))
        # TODO skip if p-value is too high?
        gts.append(snp_dosages[:, 1] + 2*snp_dosages[:, 2])
    print(f"Time: {time.time() - start}s")
    if len(gts[0]) != len(snps_to_include):
        print(var_names)
        print(snps_to_include)
        assert False

    gts = np.stack(gts, axis=1)
    assert gts.shape[1] == len(var_names)

    # rearrange the variants so that they are in sorted bp order
    # makes any subsequent visualizations useful
    sort = sorted(enumerate(var_names), key=lambda tup: int(tup[1].split('_')[1]))
    var_names = [tup[1] for tup in sort]
    sort_idxs = [tup[0] for tup in sort]
    gts = gts[:, sort_idxs]

    print("Regressing phentoypes ... ",  flush=True)
    pheno_residuals = pheno_vals - sklearn.linear_model.LinearRegression().fit(covars, pheno_vals).predict(covars)
    print("Regressing genotypes ... ",  flush=True)
    start = time.time()
    gt_residuals = np.full(gts.shape, np.nan)
    for i in range(math.ceil(gts.shape[1]/100)):
        cols = slice(i*100, (i+1)*100)
        gt_residuals[:, cols] = gts[:, cols] - sklearn.linear_model.LinearRegression().fit(covars, gts[:, cols]).predict(covars)
    print(f"Time: {time.time() - start}s")

    outdir = f'{ukb}/finemapping/susie_results/{phenotype}/{chrom}_{start_pos}_{end_pos}'
    np.save(f'{outdir}/gt_residuals.npy', gt_residuals)
    np.save(f'{outdir}/pheno_residuals.npy', pheno_residuals)
    with open(f'{outdir}/colnames.txt', 'w') as colfile:
        for var_name in var_names:
            colfile.write(var_name + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('imputation_run_name')
    parser.add_argument('phenotype')
    parser.add_argument('chrom')
    parser.add_argument('start')
    parser.add_argument('end')
    args = parser.parse_args()

    load_gts(args.imputation_run_name, args.phenotype, args.chrom, args.start, args.end)
