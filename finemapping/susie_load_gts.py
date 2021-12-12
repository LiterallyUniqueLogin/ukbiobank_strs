#!/usr/bin/env python3

import argparse
import datetime
import math
import os
import subprocess as sp
import time

import h5py
import numpy as np
import sklearn.linear_model

import load_and_filter_genotypes as lfg
import python_array_utils as utils

ukb = os.environ['UKB']

def get_str_dosages(str_dosages_dict):
    return np.sum([_len*np.sum(dosages, axis=1) for
                   _len, dosages in str_dosages_dict.items()], axis=0)

def load_gts(outreadme_fname, pheno_residuals_fname, gt_residuals_fname, outcols_fname, imputation_run_name, phenotype, chrom, start_pos, end_pos):
    filter_set_fname = f'{ukb}/finemapping/str_imp_snp_overlaps/chr{chrom}_to_filter.tab'
    outdir = f'{ukb}/finemapping/susie_results/{phenotype}/{chrom}_{start_pos}_{end_pos}'

    p_cutoff = 5e-4
   
    today = datetime.datetime.now().strftime("%Y_%M_%D")
    with open(outreadme_fname, 'w') as readme:
        readme.write(
            f'Run date: {today}\n'
            'Loading STR and SNP gts and the phenotype values, regressing out '
            'covariates from all of those, then running SuSiE. '
            f'Variants for which association tests were skipped or with p > {p_cutoff} are excluded. '
            'SNPs in the filter set are also skipped. '
            f'(Filter set at {filter_set_fname})\n'
        )

    with open(f'{ukb}/sample_qc/runs/white_brits/{phenotype}/combined_unrelated.sample') as samples_file:
        samples = np.array([line.strip() for line in samples_file][1:], dtype=int).reshape(-1, 1)

    imp_snp_samples_filepath = f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample'
    with open(imp_snp_samples_filepath) as imp_snp_samples_file:
        imp_snp_samples = np.array([line.split()[0] for line in imp_snp_samples_file][2:], dtype=int).reshape(-1, 1)

    samples_indicator = np.concatenate((samples, samples), axis=1)
    samples_merge = utils.merge_arrays(imp_snp_samples, samples_indicator)
    assert samples_merge.shape[1] == 2
    sample_idx = ~np.isnan(samples_merge[:, 1])

    pheno_covars = np.load(f'{ukb}/traits/subset_transformed_phenotypes/white_brits/{phenotype}.npy')
    shared_covars = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')
    covars = utils.merge_arrays(pheno_covars, shared_covars)

    # reorder samples to the proper order
    covars = utils.merge_arrays(samples_merge[sample_idx, 0:1], covars)
    pheno_vals = covars[:, 1]
    covars = covars[:, 2:]

    print("Regressing phentoypes ... ",  flush=True)
    pheno_residuals = pheno_vals - sklearn.linear_model.LinearRegression().fit(covars, pheno_vals).predict(covars)
    with h5py.File(pheno_residuals_fname, 'w') as pheno_residuals_file:
        pheno_residuals_dset = pheno_residuals_file.create_dataset(
            'pheno_residuals', pheno_residuals.shape, dtype='f'
        )
        pheno_residuals_dset[:] = pheno_residuals

    # first choose STRs and SNPs only with p <= p_cutoff to lessen memory burden
    print('Choosing which strs and snps to include ... ', flush=True)
    strs_to_include = set()
    snps_to_include = set()

    str_results_fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    with open(str_results_fname) as str_results:
        header = next(str_results).strip().split('\t')
    str_p_idx = header.index(f'p_{phenotype}')

    out = sp.run(
        f"grep -P '^{chrom}\t' {str_results_fname} | "
        f"awk '{{ if ({start_pos} <= $2 && $2 <= {end_pos} && ${str_p_idx + 1} <= {p_cutoff}) "
        " { print $2 } }' ",
        shell=True,
        check = True,
        capture_output=True
    )
    for line in out.stdout.decode().split('\n'):
        if line.strip() != '':
            strs_to_include.add(int(line.strip()))

    assert len(strs_to_include) != 0
    
    snps_to_filter = set()
    with open(filter_set_fname) as filter_file:
        next(filter_file) # skip header
        for line in filter_file:
            pos, ref, alt = line.strip().split('\t')[3:6]
            snps_to_filter.add((int(pos), ref, alt))

    snp_results_fname = f'{ukb}/association/results/{phenotype}/plink_snp/results.tab'
    with open(snp_results_fname) as snp_results:
        header = next(snp_results).strip().split('\t')
    snp_p_idx = header.index('P')

    out = sp.run(
        f"grep -P '^{chrom}\t' {snp_results_fname} | "
        f"awk '{{ if ({start_pos} <= $2 && $2 <= {end_pos} && ${snp_p_idx + 1} != \"NA\" && ${snp_p_idx + 1} <= {p_cutoff}) "
        " { print $2 \" \" $4 \" \" $5 } }' ",
        shell=True,
        check = True,
        capture_output=True
    )
    for line in out.stdout.decode().split('\n'):
        if line.strip() == '':
            continue
        snp = tuple(line.strip().split())
        snp = tuple([int(snp[0]), snp[1], snp[2]])
        if snp not in snps_to_filter:
            snps_to_include.add(snp)

    snp_sort_tuples = set((pos, 'SNP', ref, alt) for (pos, ref, alt) in snps_to_include)
    str_sort_tuples = set((pos, 'STR') for pos in strs_to_include)
    sorted_vars = sorted(snp_sort_tuples.union(str_sort_tuples))
    sorted_var_names =  [
        f'STR_{tuple[0]}' if tuple[1] == 'STR' else f'SNP_{tuple[0]}_{tuple[2]}_{tuple[3]}'
        for tuple in sorted_vars
    ]
    assert len(set(sorted_var_names)) == len(sorted_var_names) # make sure is unique

    print(f'# STRs: {len(strs_to_include)} # SNPs: {len(snps_to_include)}', flush=True)

    var_inclusion = { var_name : False for var_name in sorted_var_names }
    with h5py.File(gt_residuals_fname, 'w') as gt_residuals_file:
        gt_residuals_dset = gt_residuals_file.create_dataset(
            'gt_residuals',
            (pheno_covars.shape[0], len(strs_to_include) + len(snps_to_include)),
            dtype='f',
            chunks=(5, len(strs_to_include) + len(snps_to_include))
        )

        print('loading and regressing STRs... ', flush=True)
        region = f'{chrom}:{start_pos}-{end_pos}'
        str_itr = lfg.load_strs(
            imputation_run_name,
            region,
            sample_idx,
            details=False,
            var_subset=strs_to_include
        )

        start = time.time()
        for str_count, (str_dosages_dict, _, _, str_pos, str_locus_filtered, _) in enumerate(str_itr):
            str_count += 1
            print(f'loading STR {str_count} ... ', flush=True)
            str_name = f'STR_{str_pos}'
            if str_locus_filtered or str_name not in sorted_var_names:
                print(str_name, strs_to_include)
                assert False
            idx = sorted_var_names.index(str_name)
            dosages = get_str_dosages(str_dosages_dict)
            gt_residuals_dset[:, idx] = dosages - sklearn.linear_model.LinearRegression().fit(covars, dosages).predict(covars)

            assert not var_inclusion[str_name]
            var_inclusion[str_name] = True

        print(f"Time: {time.time() - start}s")
        if sum(var_inclusion.values()) != len(strs_to_include):
            print(var_inclusion)
            print(strs_to_include)
            assert False

        print('loading and regressing SNPs... ', flush=True)
        snp_itr = lfg.load_imputed_snps(
            region,
            sample_idx,
            apply_filter = False,
            details=False,
            var_subset=snps_to_include
        )
        start = time.time()
        for snp_count, (snp_dosages, alleles, _, snp_pos, snp_filtered, _) in enumerate(snp_itr):
            snp_count += 1
            if snp_count % 10 == 0:
                print(f'loading SNP {snp_count} ... ', flush=True)
            snp_name = f'SNP_{snp_pos}_{alleles[0]}_{alleles[1]}'
            if snp_filtered or snp_name not in sorted_var_names:
                print(snp_name, snps_to_include)
                assert False
            idx = sorted_var_names.index(snp_name)
            dosages = snp_dosages[:, 1] + 2*snp_dosages[:, 2]
            gt_residuals_dset[:, idx] = dosages - sklearn.linear_model.LinearRegression().fit(covars, dosages).predict(covars)
            assert not var_inclusion[snp_name]
            var_inclusion[snp_name] = True

        print(f"Time: {time.time() - start}s")
        if sum(var_inclusion.values()) != len(var_inclusion):
            print(var_inclusion)
            print(snps_to_include)
            assert False

        with open(outcols_fname, 'w') as colfile:
            for var_name in sorted_var_names:
                colfile.write(var_name + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('outreadme')
    parser.add_argument('pheno_residuals_fname')
    parser.add_argument('gt_residuals_fname')
    parser.add_argument('outcols')
    parser.add_argument('imputation_run_name')
    parser.add_argument('phenotype')
    parser.add_argument('chrom')
    parser.add_argument('start')
    parser.add_argument('end')
    args = parser.parse_args()

    load_gts(args.outreadme, args.pheno_residuals_fname, args.gt_residuals_fname, args.outcols, args.imputation_run_name, args.phenotype, args.chrom, args.start, args.end)
