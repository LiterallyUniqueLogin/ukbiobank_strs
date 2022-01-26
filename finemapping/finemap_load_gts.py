#!/usr/bin/env python3

import argparse
import os
import time

import bgen_reader
import cyvcf2
import h5py
import numpy as np

import python_array_utils as utils
import python_file_utils as file_utils
import load_and_filter_genotypes as lfg

ukb = os.environ['UKB']

def load_gts(workdir, outdir, phenotype, chrom, start_pos, end_pos, str_imputation_run_name):
    '''
    write finemap_input.master
    write gts.h5 - dataset 'gts'
    '''

    imp_snp_samples_filepath = f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample'
    with open(imp_snp_samples_filepath) as imp_snp_samples_file:
        imp_snp_samples = np.array([line.split()[0] for line in imp_snp_samples_file][2:], dtype=int)

    vcffile = f'{ukb}/str_imputed/runs/{str_imputation_run_name}/vcfs/annotated_strs/chr{chrom}.vcf.gz'
    vcf = cyvcf2.VCF(vcffile)
    str_samples = np.array([sample.split('_')[0] for sample in vcf.samples], dtype=int)
    vcf.close()

    assert np.all((imp_snp_samples == str_samples)[(str_samples > 0) & (imp_snp_samples > 0)])

    with open(f'{ukb}/sample_qc/runs/white_brits/{phenotype}/combined_unrelated.sample') as samples_file:
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
            f'{outdir}/finemap_input.z;'
            f'{outdir}/all_variants.ld;'
            f'{outdir}/finemap_output.snp;'
            f'{outdir}/finemap_output.config;'
            f'{outdir}/finemap_output.cred;'
            f'{outdir}/finemap_output.log;'
            f'{n_samples}'
        )

    with open(f'{outdir}/finemap_input.z') as finemap_input_z:
        included_strs = []
        included_imputed_snps = []
        next(finemap_input_z)
        for line in finemap_input_z:
            _id = line.split()[0]
            if _id[:4] == 'STR_':
                splits = _id.split('_')
                assert len(splits) == 2
                included_strs.append(int(splits[1]))
            elif _id[:4] == 'SNP_':
                splits = _id.split('_')
                assert len(splits) == 4
                pos, ref, alt = splits[1:]
                pos = int(pos)
                included_imputed_snps.append((pos, ref, alt))
            else:
                raise ValueError(f"ID with unknown variant type {_id}")

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

    outdir = f'{ukb}/finemapping/finemap_results/{phenotype}/{chrom}_{start_pos}_{end_pos}'

    with file_utils.temp_dir('finemap_load_gts', args) as tempdir:
        load_gts(tempdir, outdir, phenotype, chrom, start_pos, end_pos, args.str_imputation_run_name)
        file_utils.move_files(tempdir, outdir)

if __name__ == '__main__':
    main()

