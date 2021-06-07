#!/usr/bin/env python3

import argparse
import csv
import os
import pathlib
import sys
import time

import bgen_reader
import cyvcf2
import h5py
import numpy as np

import python_array_utils as utils
import load_and_filter_genotypes as lfg

ukb = os.environ['UKB']

inclusion_threshold = 0.05

def load_gts(workdir, readme, phenotype, chrom, start_pos, end_pos, str_imputation_run_name):
    '''
    write README.txt
    write finemap_input.master
    write finemap_input.z
    write gts.h5 - dataset 'gts'
    '''
    plink_results_fname = f'{ukb}/association/results/{phenotype}/plink_snp/results.tab'
    str_results_fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'

    readme.write(
        'Manually generating variant-variant LD for each imputed SNP  each STR in the region '
        'where an association was successfully '
        f'performed and had p < {inclusion_threshold}\n'
        'Correlation is STR length dosage vs SNP dosage.\n'
        'Running FINEMAP with that list of imputed SNPs and STRs.\n'
    )

    imp_snp_samples_filepath = f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample'
    with open(imp_snp_samples_filepath) as imp_snp_samples_file:
        imp_snp_samples = np.array([line.split()[0] for line in imp_snp_samples_file][2:], dtype=int)

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
    n_samples = np.sum(sample_idx)

    print(f"Working with # samples: {n_samples}", flush=True)

    with open(f'{workdir}/finemap_input.master', 'w') as finemap_master:
        finemap_master.write(
            'z;ld;snp;config;cred;log;n_samples\n'
            f'{workdir}/finemap_input.z;'
            f'{workdir}/all_variants.ld;'
            f'{workdir}/finemap_output.snp;'
            f'{workdir}/finemap_output.config;'
            f'{workdir}/finemap_output.cred;'
            f'{workdir}/finemap_output.log;'
            f'{n_samples}'
        )

    with open(f'{workdir}/finemap_input.z', 'w') as finemap_input_z:
        finemap_input_z.write('rsid chromosome position allele1 allele2 maf beta se\n')

        included_strs = []
        prev_str_pos = None
        with open(str_results_fname) as str_results_file:
            str_results_reader = csv.reader(str_results_file, delimiter='\t')
            header = next(str_results_reader)
            cols = {
                col: header.index(col) for col in
                ['chrom', 'pos', 'locus_filtered', f'p_{phenotype}', f'coeff_{phenotype}', f'se_{phenotype}']
            }

            # assumes ordered numeric chromosomes
            for result in str_results_reader:
                result_chrom = int(result[cols['chrom']])
                result_pos = int(result[cols['pos']])
                if (result_chrom, result_pos) < (chrom, start_pos):
                    continue
                if (result_chrom, result_pos) > (chrom, end_pos):
                    break
                if result[cols['locus_filtered']] != 'False' or float(result[cols[f'p_{phenotype}']]) >= inclusion_threshold:
                    continue
                if result_pos == prev_str_pos:
                    raise ValueError(f"Two STR poses at the same location {result_pos}!")
                prev_str_pos = result_pos
                included_strs.append(result_pos)

                beta = result[cols[f'coeff_{phenotype}']]
                se = result[cols[f'se_{phenotype}']]
                # I have forced there to be a unique STR per position during association testing
                # (throwing out all but the first at any single location),
                # so STR_{pos} is a unique ID
                finemap_input_z.write(
                    f'STR_{result_pos} {result_chrom:02} {result_pos} nan nan nan {beta} {se}\n'
                )
        if len(included_strs) == 0:
            pathlib.Path(f"{workdir}/finemap_output.snp").touch()
            pathlib.Path(f"{workdir}/finemap_output.config").touch()
            pathlib.Path(f"{workdir}/no_strs").touch()
            readme.write(
                "No nominally significant (p<=0.05) STRs were found in the region, "
                "so finemapping is being skipped."
            )
            print(
                "No nominally significant (p<=0.05) STRs were found in the region, "
                "so finemapping is being skipped.",
                flush = True
            )
            sys.exit()

        included_imputed_snps = []
        with open(plink_results_fname) as plink_result_file:
            plink_results_reader = csv.reader(plink_result_file, delimiter='\t')
            header = next(plink_results_reader)
            cols = {
                col: header.index(col) for col in
                ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'BETA', 'SE', 'P', 'ERRCODE']
            }

            for result in plink_results_reader:
                result_chrom = int(result[cols['#CHROM']])
                result_pos = int(result[cols['POS']])
                ref = result[cols['REF']]
                alt = result[cols['ALT']]
                if (result_chrom, result_pos) < (chrom, start_pos):
                    continue
                if (result_chrom, result_pos) > (chrom, end_pos):
                    break
                if result[cols['ERRCODE']] != '.' or float(result[cols['P']]) >= inclusion_threshold:
                    continue
                # snps can only be uniquely identified by pos, ref, alt
                # some IDs are duplicate at the same or even different locations
                # so those don't help
                included_imputed_snps.append((result_pos, ref, alt))

                beta = result[cols['BETA']]
                se = result[cols['SE']]
                finemap_input_z.write(
                    f'SNP_{result_pos}_{ref}_{alt} {result_chrom:02} {result_pos} {ref} {alt} nan {beta} {se}\n'
                )

    n_snps = len(included_imputed_snps)
    n_strs = len(included_strs)

    n_variants = n_snps + n_strs

    print(f'Num snps: {n_snps}, num strs: {n_strs}, total: {n_variants}', flush=True)

    chunk_len = 2**6
    with h5py.File(f'{workdir}/gts.h5', 'w') as h5file:
        if n_variants >= chunk_len:
            gts = h5file.create_dataset("gts", (n_variants, n_samples), chunks=(chunk_len, n_samples))
        else:
            gts = h5file.create_dataset("gts", (n_variants, n_samples))

        print("Loading strs ...", flush=True)
        start = time.time()
        next_strs = iter(included_strs)
        next_str = next(next_strs)
        region = f'{chrom}:{start_pos}-{end_pos}'
        all_str_itr = lfg.load_strs(str_imputation_run_name, region, sample_idx)
        next(all_str_itr) # skip locus_details list
        str_num = 0
        for str_dosages_dict, _, _, str_pos, str_locus_filtered, _ in all_str_itr:
            if str_locus_filtered:
                continue
            if str_pos < next_str:
                continue
            elif str_pos > next_str:
                raise ValueError(f"Skipped str {next_str}!")
            str_dosages = get_str_dosages(str_dosages_dict)

            gts[str_num, :] = str_dosages

            str_num += 1
            if str_num % 5 == 0:
                print(f'{str_num}/{n_strs} loaded. ETA: {(time.time() - start)*(n_strs-str_num)/str_num}sec', flush=True)

            try:
                next_str = next(next_strs)
            except StopIteration:
                break
        print(f'Done loading STRs. Time: {time.time() - start}sec', flush=True)

        print("Loading snps ...", flush=True)
        snp_bgen = bgen_reader.open_bgen(
            f'{ukb}/array_imputed/ukb_imp_chr{chrom}_v3.bgen',
            verbose=False
        )

        start = time.time()
        next_snps = iter(included_imputed_snps)
        next_snp_pos, next_snp_ref, next_snp_alt = next(next_snps)
        snp_num = 0
        for all_snp_num, snp_pos in enumerate(snp_bgen.positions):
            if snp_pos < next_snp_pos:
                continue
            elif snp_pos > next_snp_pos:
                raise ValueError(f"Skipped snp pos {next_snp_pos}!")

            # match not just on pos but also ref, alt
            # assume that multiple SNPs at the same pos have the same ordering
            ref, alt = snp_bgen.allele_ids[all_snp_num].split(',')
            if (ref, alt) != (next_snp_ref, next_snp_alt):
                continue

            snp_probs = snp_bgen.read(all_snp_num).squeeze()[sample_idx, :]
            assert len(snp_probs.shape) == 2
            assert snp_probs.shape[1] == 3
            snp_dosages = snp_probs[:, 1] + 2*snp_probs[:, 2]
            gts[n_strs + snp_num, :] = snp_dosages

            snp_num += 1
            if snp_num % 100 == 0:
                print(f'{snp_num}/{n_snps} loaded. ETA: {(time.time() - start)*(n_snps-snp_num)/snp_num}sec', flush=True)

            try:
                next_snp_pos, next_snp_ref, next_snp_alt = next(next_snps)
            except StopIteration:
                break
        print(f'Done loading SNPs. Time: {time.time() - start}sec', flush=True)

def get_str_dosages(str_dosages_dict):
    return np.sum([_len*np.sum(dosages, axis=1) for
                   _len, dosages in str_dosages_dict.items()], axis=0)

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
        load_gts(workdir, readme, phenotype, chrom, start_pos, end_pos, args.str_imputation_run_name)

if __name__ == '__main__':
    main()

