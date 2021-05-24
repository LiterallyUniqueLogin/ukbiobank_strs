#!/usr/bin/env python3

import argparse
import csv
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

inclusion_threshold = 0.005

def write_ldstore_input(readme, workdir, phenotype, chrom, start_pos, end_pos):
    plink_results_fname = f'{ukb}/association/results/{phenotype}/plink_snp/results.tab'
    readme.write(
        'Running ldstore2 to get pairwise LD for each imputed SNP in the region '
        'that plink ran an association for without error and that reached p < '
        f'{inclusion_threshold}\n'
    )
    included_imputed_snps = []
    with open(plink_results_fname) as plink_result_file:
        plink_results_reader = csv.reader(plink_result_file, delimiter='\t')
        header = next(plink_results_reader)
        cols = {
            col: header.index(col) for col in
            ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'BETA', 'SE', 'P', 'ERRCODE']
        }

        with open(f'{workdir}/ldstore_input.z', 'w') as ldstore_input, \
                open(f'{workdir}/finemap_input.z', 'w') as finemap_input_z:
            ldstore_input.write('rsid chromosome position allele1 allele2\n')
            finemap_input_z.write('rsid chromosome position allele1 allele2 beta se\n')
            # assumes ordered numeric chromosomes
            for result in plink_results_reader:
                result_chrom = int(result[cols['#CHROM']])
                result_pos = int(result[cols['POS']])
                if (result_chrom, result_pos) < (chrom, start_pos):
                    continue
                if (result_chrom, result_pos) > (chrom, end_pos):
                    break
                if result[cols['ERRCODE']] != '.' or float(result[cols['P']]) >= inclusion_threshold:
                    continue
                included_imputed_snps.append(result_pos)

                rsid = result[cols['ID']]
                allele1 = result[cols['REF']]
                allele2 = result[cols['ALT']]
                beta = result[cols['BETA']]
                se = result[cols['SE']]
                ldstore_input.write(
                    f'{rsid} {result_chrom:02} {result_pos} {allele1} {allele2}\n'
                )
                finemap_input_z.write(
                    f'{rsid} {result_chrom:02} {result_pos} nan nan {beta} {se}\n'
                )

    outsample_fname = f'{workdir}/ldstore_input.incl'
    nsamples = 0
    with open(f'{ukb}/sample_qc/runs/{phenotype}/combined_unrelated.sample') as sample_file, \
            open(outsample_fname, 'w') as outsample_file:
        first = True
        for line in sample_file:
            if first:
                first = False
                continue
            nsamples += 1
            outsample_file.write(line)

    bgen_fname = f'{ukb}/array_imputed/ukb_imp_chr{chrom}_v3.bgen'
    with open(f'{workdir}/ldstore_input.master', 'w') as ldstore_master:
        ldstore_master.write('z;bgen;bgi;sample;incl;bdose;ld;n_samples\n')
        ldstore_master.write(
            f'{workdir}/ldstore_input.z;'
            f'{bgen_fname};'
            f'{bgen_fname}.bgi;'
            f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample;'
            f'{outsample_fname};'
            f'{workdir}/ldstore_temp.bdose;'
            f'{workdir}/ldstore_out.ld;'
            f'{nsamples}'
        )

    return included_imputed_snps

def run_ldstore(workdir):
    out = sp.run(
        f'{ukb}/utilities/ldstore/ldstore_v2.0_x86_64 --write-text '
        f'--in-files {workdir}/ldstore_input.master '
        f'--write-bdose --bdose-version 1.0 '
        '--memory 19.9 --n-threads 8',
        shell=True,
        capture_output=True
    )
    print(out.stdout.decode(), flush=True)
    print(out.stderr.decode(), file=sys.stderr, flush=True)
    out.check_returncode()

def prep_finemap_input(
        readme,
        workdir,
        phenotype,
        chrom,
        start_pos,
        end_pos,
        str_imputation_run_name,
        included_imputed_snps):

    snp_ld_matrix = np.genfromtxt(f'{workdir}/ldstore_out.ld')
    assert len(snp_ld_matrix.shape) == 2
    assert snp_ld_matrix.shape[0] == snp_ld_matrix.shape[1]
    assert snp_ld_matrix.shape[0] == len(included_imputed_snps)
    num_snps = snp_ld_matrix.shape[0]

    readme.write(
        'Manually generating STR-STR LD and STR-SNP LD for each imputed SNP included in '
        'the ldstore2 run and each STR in the region where an association was successfully '
        f'performed and had p < {inclusion_threshold}\n'
        'Running FINEMAP with that list of imputed SNPs and STRs.\n'
    )

    str_results_fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    included_strs = []
    with open(str_results_fname) as str_result_file, \
            open(f'{workdir}/finemap_input.z', 'a') as finemap_input_z:
        str_results_reader = csv.reader(str_result_file, delimiter='\t')
        header = next(str_results_reader)
        cols = {
            col: header.index(col) for col in
            ['chrom', 'pos', 'locus_filtered', f'p_{phenotype}', f'coeff_{phenotype}', f'se_{phenotype}']
        }

        last_str_pos = None
        num_strs = 0
        for result in str_results_reader:
            num_strs += 1
            result_chrom = int(result[cols['chrom']])
            result_pos = int(result[cols['pos']])
            if (result_chrom, result_pos) < (chrom, start_pos):
                continue
            if (result_chrom, result_pos) > (chrom, end_pos):
                break
            if result[cols['locus_filtered']] != 'False' or float(result[cols[f'p_{phenotype}']]) >= inclusion_threshold:
                continue
            included_strs.append(result_pos)

            if result_pos == last_str_pos:
                raise ValueError(f"Already encountered pos {result_pos}!")

            beta = result[cols[f'coeff_{phenotype}']]
            se = result[cols[f'coeff_{phenotype}']]
            finemap_input_z.write(
                f'STR_{result_pos} {result_chrom} {result_pos} nan nan nan {beta} {se}\n'
            )

    imp_snp_samples_filepath = f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample'
    with open(imp_snp_samples_filepath) as imp_snp_samples_file:
        imp_snp_samples = np.array([line.split()[0] for line in imp_snp_samples_file][2:], dtype=int)

    region = f'{chrom}:{start_pos}-{end_pos}'

    vcffile = f'{ukb}/str_imputed/runs/{str_imputation_run_name}/vcfs/annotated_strs/chr{chrom}.vcf.gz'
    vcf = cyvcf2.VCF(vcffile)
    str_samples = np.array([sample.split('_')[0] for sample in vcf.samples], dtype=int)
    vcf.close()

    assert np.all((imp_snp_samples == str_samples)[(str_samples > 0) & (imp_snp_samples > 0)])

    with open(f'{ukb}/sample_qc/runs/{phenotype}/combined_unrelated.sample') as samples_file:
        samples = np.array([line.strip() for line in samples_file][1:], dtype=int).reshape(-1, 1)

    samples_indicator = np.concatenate((samples, samples), axis=1)
    samples_merge = utils.merge_arrays(str_samples.reshape(-1, 1), samples_indicator)
    assert samples_merge.shape[1] == 2
    sample_idx = ~np.isnan(samples_merge[:, 1])

    num_variants = num_strs + num_snps
    print(f'Num snps: {num_snps}, num strs: {num_strs}, total: {num_variants}', flush=True)
    total_corrs = num_strs*num_snps + num_strs**2

    full_ld_matrix = np.full((num_variants, num_variants), np.nan)
    full_ld_matrix[:num_snps, :num_snps] = snp_ld_matrix

    snp_bgen = bgen_reader.open_bgen(
        f'{ukb}/array_imputed/ukb_imp_chr{chrom}_v3.bgen',
        verbose=False
    )

    corrs_complete = 0
    start = time.time()
    next_strs = enumerate(included_strs)
    str_num, next_str = next(next_strs)
    all_str_itr = lfg.load_strs(str_imputation_run_name, region, sample_idx)
    next(all_str_itr) # skip locus_details list
    done_done = False
    for str_dosages_dict, _, _, str_pos, str_locus_filtered, _ in all_str_itr:
        if str_locus_filtered:
            continue
        if str_pos < next_str:
            continue
        elif str_pos == next_str:
            try:
                str_num, next_str = next(next_strs)
            except StopIteration:
                done_done = True
        else:
            raise ValueError(f"Skipped str {next_str}!")
        str_dosages = get_str_dosages(str_dosages_dict)

        next_snps = enumerate(included_imputed_snps)
        snp_num, next_snp = next(next_snps)
        snp_done = False
        for all_snp_num, snp_pos in enumerate(snp_bgen.positions):
            if snp_pos < next_snp:
                continue
            elif snp_pos == next_snp:
                try:
                    snp_num, next_snp = next(next_snps)
                except StopIteration:
                    snp_done = True
            else:
                raise ValueError(f"Skipped snp {next_snp}!")
            snp_probs = snp_bgen.read(all_snp_num).squeeze()[sample_idx, :]
            assert len(snp_probs.shape) == 2
            assert snp_probs.shape[1] == 3
            snp_dosages = snp_probs[:, 1] + 2*snp_probs[:, 2]
            corr = np.corrcoef(str_dosages, snp_dosages)[0,1]
            full_ld_matrix[snp_num, num_snps + str_num] = corr
            full_ld_matrix[num_snps + str_num, snp_num] = corr
            corrs_complete += 1
            if corrs_complete % 20 == 0:
                print(f'{corrs_complete}/{total_corrs} correlation. ETA: {(time.time() - start)*(total_corrs-corrs_complete)/corrs_complete}sec', flush=True)
            if snp_done:
                break

        next_strs2 = enumerate(included_strs)
        str_additional_num, next_str2 = next(next_strs2)
        all_str_itr2 = lfg.load_strs(str_imputation_run_name, f'{chrom}:{str_pos}-{end_pos}', sample_idx)
        next(all_str_itr2) # skip locus_details list
        str_done = False
        for str_dosages_dict2, _, _, str_pos2, str_locus_filtered2, _ in all_str_itr2:
            if str_locus_filtered2:
                continue
            if str_pos2 < next_str2:
                continue
            elif str_pos2 == next_str2:
                try:
                    str_additional_num, next_str2 = next(next_strs2)
                except StopIteration:
                    str_done = True
            else:
                raise ValueError(f"Skipped str {next_str2}!")
            str_dosages2 = get_str_dosages(str_dosages_dict2)
            corr = np.corrcoef(str_dosages, str_dosages2)[0,1]
            full_ld_matrix[num_snps + str_num, num_snps + str_num + str_additional_num] = corr
            full_ld_matrix[num_snps + str_num + str_additional_num, num_snps + str_num] = corr
            corrs_complete += 1
            if corrs_complete % 20 == 0:
                print(f'{corrs_complete}/{total_corrs} correlation. ETA: {(time.time() - start)*(total_corrs-corrs_complete)/corrs_complete}sec', flush=True)
            if str_done:
                break
        if done_done:
            break

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
        '--n-configs-top 100 '
        '--n-threads 8',
        shell=True,
        capture_output=True
    )

    print(out.stdout.decode(), flush=True)
    print(out.stderr.decode(), file=sys.stderr, flush=True)
    out.check_returncode()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('chrom', type=int)
    parser.add_argument('start_pos', type=int)
    parser.add_argument('end_pos', type=int)
    parser.add_argument('str_imputation_run_name')
    args = parser.parse_args()

    phenotype = args.phenotype
    chrom = args.chrom
    start_pos = args.start_pos
    end_pos = args.end_pos
    assert start_pos < end_pos

    workdir = f'{ukb}/finemapping/finemap_results/{phenotype}/{chrom}_{start_pos}_{end_pos}'

    with open(f'{workdir}/README.txt', 'w') as readme:
        start = time.time()
        print('Writing ldstore input ...', flush=True)
        included_imputed_snps = write_ldstore_input(readme, workdir, phenotype, chrom, start_pos, end_pos)
        print(f'Done writing ldstore input. Elapsed time = {time.time() - start}sec', flush=True)

        start = time.time()
        print('Running ldstore ... ', flush=True)
        run_ldstore(workdir)
        print(f'Done running ldstore. Elapsed time = {time.time() - start}sec', flush=True)

        start = time.time()
        print('Prepping FINEMAP input ... ', flush=True)
        prep_finemap_input(
            readme,
            workdir,
            phenotype,
            chrom,
            start_pos,
            end_pos,
            args.str_imputation_run_name,
            included_imputed_snps
        )
        print(f'Done prepping FINEMAP input. Elapsed time = {time.time() - start}sec', flush=True)

        start = time.time()
        print('Running FINEMAP ... ', flush=True)
        run_finemap(workdir)
        print(f'Done running FINEMAP. Elapsed time = {time.time() - start}sec', flush=True)

if __name__ == '__main__':
    main()
