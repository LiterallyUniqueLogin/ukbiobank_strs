#!/usr/bin/env python3

import argparse
import os
import subprocess as sp
import sys
import time

import bgen_reader
import cyvcf2
import numpy as np
import pandas as pd

import python_array_utils as utils
import load_and_filter_genotypes as lfg

ukb = os.environ['UKB']

def write_ldstore_input(workdir, phenotype, chrom, start, end):
    plink_results_fname = f'{ukb}/association/results/{phenotype}/plink_snp/results.tab'
    plink_results = utils.df_to_recarray(pd.read_csv(
        plink_results_fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        usecols=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'BETA', 'SE', 'ERRCODE'],
        dtype={
            '#CHROM': int,
            'POS': int,
            'ID': object,
            'REF': object,
            'ALT': object,
            'BETA': float,
            'SE': float,
            'ERRCODE': object
        }
    ))

    snp_rsids = []
    rsid_set = set()

    with open(f'{workdir}/ldstore_input.z', 'w') as ldstore_input, \
            open(f'{workdir}/finemap_input.z', 'w') as finemap_input_z:
        ldstore_input.write('rsid chromosome position allele1 allele2\n')
        # assumes ordered numeric chromosomes
        for idx in range(plink_results.shape[0]):
            result = plink_results[idx]
            result_chrom = result['#CHROM']
            result_pos = result['POS']
            if (result_chrom, result_pos) < (chrom, start):
                continue
            if (result_chrom, result_pos) > (chrom, end):
                break
            if result['ERRCODE'] != '.':
                continue

            rsid = result['ID']
            if rsid in rsid_set:
                raise ValueError(f"Already encountered rsid {rsid}!")
            rsid_set.add(rsid)
            snp_rsids.append(rsid)
            allele1 = result['REF']
            allele2 = result['ALT']
            beta = result['BETA']
            se = result['SE']
            ldstore_input.write(
                f'{rsid} {result_chrom} {result_pos} {allele1} {allele2}\n'
            )
            finemap_input_z.write(
                f'{rsid} {result_chrom} {result_pos} nan nan nan {beta} {se}\n'
            )

    outsample_fname = f'{workdir}/ldstore_input.sample'
    with open(f'{ukb}/sample_qc/runs/{phenotype}/combined_unrelated.sample') as sample_file, \
            open(outsample_fname, 'w') as outsample_file:
        first = True
        for line in sample_file:
            if first:
                first = False
                continue
            outsample_file.write(line)

    bgen_fname = f'{ukb}/array_imputed/ukb_imp_chr{chrom}_v3.bgen'
    with open(f'{workdir}/ldstore_input.master', 'w') as ldstore_master:
        ldstore_master.write('z;bgen;bgi;sample;incl;bdose;ld\n')
        ldstore_master.write(
            f'{workdir}/ldstore_input.z;'
            f'{bgen_fname};'
            f'{bgen_fname}.bgi;'
            f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample;'
            f'{outsample_fname};'
            f'{workdir}/ldstore_temp.bdose;',
            f'{workdir}/ldstore_out.ld'
        )

    return snp_rsids

def run_ldstore(workdir):
    out = sp.run(
        f'{ukb}/utilities/ldstore/ldstore_v2.0_x86_64 --write-text '
        f'--in-files {workdir}/ldstore_input.master '
        f'--write-bdose --bdose-version 1.0 '
        '--memory 1.9Gb ',
        shell=True,
        capture_output=True
    )

    if out.returncode != 0:
        print(out.stdout.decode())
        print(out.stderr.decode(), file=sys.stderr)
        out.check_returncode()

def prep_finemap_input(
        workdir,
        phenotype,
        chrom,
        start,
        end,
        str_imputation_run_name,
        snp_rsids):

    snp_ld_matrix = np.genfromtxt(f'{workdir}/ldstore_out.ld')
    assert snp_ld_matrix.shape == 2
    assert snp_ld_matrix.shape[0] == snp_ld_matrix.shape[1]
    num_snps = snp_ld_matrix.shape[0]

    str_results_fname = f'{ukb}/association/results/{phenotype}/plink_snp/results.tab'
    str_results = utils.df_to_recarray(pd.read_csv(
        str_results_fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        usecols=['chrom', 'pos', 'locus_filtered', f'coeff_{phenotype}', f'se_{phenotype}'],
        dtype={
            'chrom': int,
            'pos': int,
            'locus_filtered': object,
            f'coeff_{phenotype}': float,
            f'se_{phenotype}': float
        }
    ))

    str_poses = []

    with open(f'{workdir}/finemap_input.z', 'a') as finemap_input_z:
        for idx in range(str_results.shape[0]):
            result = str_results[idx]
            result_chrom = result['chrom']
            result_pos = result['pos']
            if (result_chrom, result_pos) < (chrom, start):
                continue
            if (result_chrom, result_pos) > (chrom, end):
                break
            if result['locus_filtered'] != 'False':
                continue

            if str_poses[-1] == result_pos:
                raise ValueError(f"Already encountered pos {result_pos}!")
            str_poses.append(result_pos)

            beta = result[f'coeff_{phenotype}']
            se = result[f'coeff_{phenotype}']
            finemap_input_z.write(
                f'STR_{result_pos} {result_chrom} {result_pos} nan nan nan {beta} {se}\n'
            )
    num_strs = len(str_poses)

    samples_filepath = f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample;'
    with open(samples_filepath) as all_samples_file:
        all_samples = [line.split()[0] for line in all_samples_file][2:]

    region = f'{chrom}:{start}-{end}'

    vcffile = f'{ukb}/str_imputed/runs/{str_imputation_run_name}/vcfs/annotated_strs/chr{chrom}.vcf.gz'
    vcf = cyvcf2.VCF(vcffile)
    assert [sample.split('_')[0] for sample in vcf.samples] == all_samples
    vcf.close()

    with open(f'{ukb}/sample_qc/runs/{phenotype}/combined_unrelated.sample') as samples_file:
        samples = np.array([line.strip() for line in samples_file][1:]).reshape(-1, 1)

    samples_indicator = np.concatenate((samples, samples), axis=1)
    samples_merge = utils.merge_arrays(all_samples, samples_indicator)
    assert samples_merge.shape[1] == 2
    sample_idx = ~np.isnan(samples_merge[:, 1])

    num_variants = num_strs + num_snps
    print(f'Num snps: {num_snps}, num strs: {num_strs} total: {num_variants}', flush=True)
    total_corrs = num_strs*num_snps + num_strs**2

    full_ld_matrix = np.full((num_variants, num_variants), np.nan)
    full_ld_matrix[:num_snps, :num_snps] = snp_ld_matrix

    corrs_complete = 0
    start = time.time()
    for str_num, (str_dosages_dict, _, _, str_pos, str_locus_filtered, _) in enumerate(lfg.load_strs(
            str_imputation_run_name,
            region,
            sample_idx)):
        if str_locus_filtered:
            continue
        str_dosages = get_str_dosages(str_dosages_dict)

        for snp_num, (snp_dosage_array, _, _, _, snp_locus_filtered, _) in \
                enumerate(lfg.load_imputed_snps(region, sample_idx)):
            if snp_locus_filtered:
                continue
            snp_dosages = snp_dosage_array[:, 1] + 2*snp_dosage_array[:, 2]
            corr = np.corrcoef(str_dosages, snp_dosages)[0,1]
            full_ld_matrix[snp_num, num_snps + str_num] = corr
            full_ld_matrix[num_snps + str_num, snp_num] = corr
            corrs_complete += 1
            if corrs_complete % 20 == 0:
                print(f'{corrs_complete}/{total_corrs} correlation. ETA: {(time.time() - start)*(total_corrs-corrs_complete)/corrs_complete}sec', flush=True)

        for str_additional_num, (str_dosages_dict2, _, _, _, str_locus_filtered2, _) in enumerate(lfg.load_strs(
                str_imputation_run_name,
                f'{chrom}:{str_pos}-{end}',
                sample_idx)):
            if str_locus_filtered2:
                continue
            str_dosages2 = get_str_dosages(str_dosages_dict2)
            corr = np.corrcoef(str_dosages, str_dosages2)[0,1]
            full_ld_matrix[num_snps + str_num, num_snps + str_num + str_additional_num] = corr
            full_ld_matrix[num_snps + str_num + str_additional_num, num_snps + str_num] = corr
            corrs_complete += 1
            if corrs_complete % 20 == 0:
                print(f'{corrs_complete}/{total_corrs} correlation. ETA: {(time.time() - start)*(total_corrs-corrs_complete)/corrs_complete}sec', flush=True)

    np.savetxt(f'{workdir}/all_variants.ld', full_ld_matrix)

    with open(f'{workdir}/finemap_input.master', 'w') as finemap_master:
        finemap_master.write(
            'z;ld;snp;config;cred;log;n_samples\n'
            f'{workdir}/finemap_input.z;'
            f'{workdir}/all_variants.ld;'
            f'{workdir}/finemap_output.snp;'
            f'{workdir}/finemap_output.config;'
            f'{workdir}/finemap_output.cred;'
            f'{workdir}/finemap_output.log'
        )


def get_str_dosages(str_dosages_dict):
    return np.sum([_len*np.sum(dosages, axis=1) for
                   _len, dosages in str_dosages_dict.items()], axis=0)

def run_finemap(workdir):
    out = sp.run(
        f'{ukb}/utilities/finemap/finemap_v1.4_x86_64 --sss '
        f'--in-files {workdir}/finemap_input.master '
        '--log '
        '--n-configs-top 100 ',
        shell=True,
        capture_output=True
    )

    if out.returncode != 0:
        print(out.stdout.decode())
        print(out.stderr.decode(), file=sys.stderr)
        out.check_returncode()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('chrom')
    parser.add_argument('start', type=int)
    parser.add_argument('end', type=int)
    parser.add_argument('str-imputation-run_name')
    args = parser.parse_args()

    phenotype = args.phenotype
    chrom = args.chrom
    start = args.start
    end = args.end
    assert start < end

    workdir = f'{ukb}/finemapping/finemap_results/{phenotype}/{chrom}_{start}_{end}'

    start = time.time()
    print('Writing ldstore input ...', flush=True)
    snp_rsids = write_ldstore_input(workdir, phenotype, chrom, start, end)
    print(f'done. Elapsed time = {time.time() - start}sec', flush=True)

    start = time.time()
    print('Running ldstore ... ', flush=True)
    run_ldstore(workdir)
    print(f'Done running ldstore. Elapsed time = {time.time() - start}sec', flush=True)

    start = time.time()
    print('Prepping finemap input ... ', flush=True)
    prep_finemap_input(
        workdir,
        phenotype,
        chrom,
        start,
        end,
        args.str_imputation_run_name,
        snp_rsids
    )
    print(f'Done prepping finemap input. Elapsed time = {time.time() - start}sec', flush=True)

if __name__ == '__main__':
    main()
