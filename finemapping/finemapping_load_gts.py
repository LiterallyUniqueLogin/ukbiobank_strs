#!/usr/bin/env python3

import argparse
import datetime
import shutil
import tempfile
import time

import h5py
import numpy as np
import polars as pl
import sklearn.linear_model

import load_and_filter_genotypes as lfg
import python_array_utils as utils
import sample_utils
#import load_PACSIN2

def get_str_dosages(str_dosages_dict):
    return np.sum([_len*np.sum(dosages, axis=1) for
                   _len, dosages in str_dosages_dict.items()], axis=0)

def store_values(gts_dset, stored_indexes, gt_temp, covars, take_residuals):
    assert not np.any(np.isnan(gt_temp))
    if take_residuals:
        residuals = gt_temp - sklearn.linear_model.LinearRegression().fit(covars, gt_temp).predict(covars)
        for count, idx in enumerate(stored_indexes):
            gts_dset[:, idx] = residuals[:, count]
    else:
        for count, idx in enumerate(stored_indexes):
            gts_dset[:, idx] = gt_temp[:, count]

def load_gts(readme_fname, gts_fname, str_vcf, snp_bgen, varname_fname, all_samples_fname, phenotype_samples_fname, phenotype, chrom, start_pos, end_pos, varname_header, pheno_fname, shared_covars_fname, pheno_out, hardcalls):
    today = datetime.datetime.now().strftime("%Y_%M_%D")
    with open(readme_fname, 'w') as readme:
        readme.write(f'Run date: {today}\nLoading STR and SNP gts from {varname_fname}.\n')
        if pheno_fname:
            readme.write('Regressing covariates out.\n')
            readme.write('Regressing covariates out of phenotype as well.\n')
        if not hardcalls:
            readme.write('Using dosages.\n')
        else:
            readme.write('Using hardcalls.\n')

    print('Loading samples (and phenotypes if selected) ... ', flush=True)
    sample_idx = sample_utils.get_samples_idx(all_samples_fname, phenotype_samples_fname)
    n_samples = np.sum(sample_idx)

    if pheno_fname:
        pheno_covars = np.load(pheno_fname)#f'{ukb}/traits/subset_transformed_phenotypes/white_brits/{phenotype}.npy')
        shared_covars = np.load(shared_covars_fname)#f'{ukb}/traits/shared_covars/shared_covars.npy')
        covars = utils.merge_arrays(pheno_covars, shared_covars)

        # reorder samples to the proper order
        ordered_samples = sample_utils.get_ordered_samples(all_samples_fname, phenotype_samples_fname)
        covars = utils.merge_arrays(ordered_samples.reshape(-1, 1), covars)
        pheno_vals = covars[:, 1]
        covars = covars[:, 2:]

        print("Regressing phentoypes ... ",  flush=True)
        pheno_residuals = pheno_vals - sklearn.linear_model.LinearRegression().fit(covars, pheno_vals).predict(covars)
        with h5py.File(pheno_out, 'w') as pheno_residuals_file:
            pheno_residuals_dset = pheno_residuals_file.create_dataset(
                'pheno_residuals', pheno_residuals.shape, dtype='f'
            )
            pheno_residuals_dset[:] = pheno_residuals
    else:
        covars = None

    print('Loading the strs and snps to include ... ', flush=True)
    varnames = pl.scan_csv(
        varname_fname,
        sep=' ',
        has_header = varname_header,
        with_column_names = lambda colnames: ['varnames'] + colnames[1:]
    ).collect()['varnames'].to_list()
    
    strs_to_include = set()
    snps_to_include = set()
    pacsin2_vars = set()
    for varname in varnames:
        if varname[:3] == 'STR':
            strs_to_include.add(int(varname.split('_', 1)[1]))
        elif varname[:3] == 'SNP':
            splits = varname.split('_')
            assert len(splits) == 4
            snps_to_include.add((int(splits[1]), splits[2], splits[3]))
        else:
            assert varname[:7] == 'PACSIN2'
            pacsin2_vars.add(int(varname.split('_', 2)[2]))
    print(f'# STRs: {len(strs_to_include)} # SNPs: {len(snps_to_include)}', flush=True)

    var_inclusion = { var_name : False for var_name in varnames }

    region = f'{chrom}:{start_pos}-{end_pos}'
    chunk_len = 2**6
    with h5py.File(gts_fname, 'w') as gts_file:
        if len(varnames) >= chunk_len:
            gts_dset = gts_file.create_dataset(
                'gts' if not pheno_fname else 'gt_residuals',
                (n_samples, len(varnames)),
                dtype='f',
                chunks=(n_samples, chunk_len)
            )
        else:
            gts_dset = gts_file.create_dataset(
                'gts' if not pheno_fname else 'gt_residuals',
                (n_samples, len(varnames)),
                dtype='f',
            )

        n_temp = 10
        gt_temp = np.full((n_samples, n_temp), np.nan)
        stored_indexes = []

        if len(pacsin2_vars) > 0:
            if not pheno_fname:
                print('loading PACSIN2 STRs ... ', flush=True)
            else:
                print('loading and regressing PACSIN2 STRs ... ', flush=True)
            PACSIN2_count = 0
            PACSIN2_itr = load_PACSIN2.get_gt_itr(sample_idx)
            next(PACSIN2_itr) # skip details
            for (PACSIN2_dosages_or_hardcalls, _, _, PACSIN2_pos, _, _) in PACSIN2_itr:
                if PACSIN2_pos not in pacsin2_vars:
                    continue
                PACSIN2_count += 1
                PACSIN2_name = f'PACSIN2_STR_{PACSIN2_pos}'
                assert PACSIN2_name in varnames
                del var_inclusion[PACSIN2_name]
                idx = varnames.index(PACSIN2_name)
                stored_indexes.append(idx)
                if not hardcalls:
                    gts = get_str_dosages(PACSIN2_dosages_or_hardcalls)
                else:
                    gts = np.sum(PACSIN2_dosages_or_hardcalls, axis=1)

                gt_temp[:, (PACSIN2_count - 1) % n_temp] = gts

            assert np.all(np.isnan(gt_temp[:, len(stored_indexes):]))
            gt_temp = gt_temp[:, :len(stored_indexes)]
            store_values(gts_dset, stored_indexes, gt_temp, covars, pheno_fname)
            print('done!', flush=True)
            gt_temp = np.full((n_samples, n_temp), np.nan)
            stored_indexes = []

        if pheno_fname:
            print('loading and regressing STRs... ', flush=True)
        else:
            print('loading STRs... ', flush=True)

        start = time.time()
        str_itr = lfg.load_strs(
            str_vcf,
            region,
            sample_idx,
            details=False,
            var_subset=strs_to_include,
            hardcalls=hardcalls,
            both_poses = True
        )
        for str_count, (str_dosages_or_hardcalls, _, _, (str_start_pos, str_written_pos), str_locus_filtered, _) in enumerate(str_itr):
            str_count += 1
            print(f'loading STR {str_count}, time/STR: {(time.time() - start)/str_count:.2}s ... ', flush=True)
            done = False
            for pos in (str_start_pos, str_written_pos):
                str_name = f'STR_{pos}'
                if str_name in varnames:
                    done = True
                    break
            if not done or str_locus_filtered or var_inclusion[str_name]:
                print(str_start_pos, str_written_pos, varnames, var_inclusion[str_name])
                assert False

            var_inclusion[str_name] = True

            idx = varnames.index(str_name)
            stored_indexes.append(idx)
            if not hardcalls:
                gts = get_str_dosages(str_dosages_or_hardcalls)
            else:
                gts = np.sum(str_dosages_or_hardcalls, axis=1)

            gt_temp[:, (str_count - 1) % n_temp] = gts
            if str_count % n_temp == 0:
                assert len(stored_indexes) == n_temp
                store_values(gts_dset, stored_indexes, gt_temp, covars, pheno_fname)
                gt_temp = np.full((n_samples, n_temp), np.nan)
                stored_indexes = []

        if len(stored_indexes) > 0:
            assert np.all(np.isnan(gt_temp[:, len(stored_indexes):]))
            gt_temp = gt_temp[:, :len(stored_indexes)]
            store_values(gts_dset, stored_indexes, gt_temp, covars, pheno_fname)
        gt_temp = np.full((n_samples, n_temp), np.nan)
        stored_indexes = []

        print(f"Time: {time.time() - start}s")
        if sum(var_inclusion.values()) != len(strs_to_include):
            print(var_inclusion)
            print(strs_to_include)
            assert False

        if not pheno_fname:
            print('loading SNPs... ', flush=True)
        else:
            print('loading and regressing SNPs... ', flush=True)

        snp_itr = lfg.load_imputed_snps(
            snp_bgen,
            None,
            region,
            sample_idx,
            apply_filter=False,
            details=False,
            var_subset=snps_to_include,
            hardcalls = hardcalls
        )
        for snp_count, (snp_dosages_or_hardcalls, alleles, _, snp_pos, snp_filtered, _) in enumerate(snp_itr):
            snp_count += 1
            if snp_count % n_temp == 0:
                if snp_count > 3:
                    print(f'loading SNP {snp_count}, time/snp: {(time.time() - start)/snp_count:.2}s ... ', flush=True)
                else:
                    print(f'loading SNP {snp_count} ... ', flush=True)
            if snp_count == 3:
                start = time.time()
            snp_name = f'SNP_{snp_pos}_{alleles[0]}_{alleles[1]}'
            if snp_filtered or snp_name not in varnames:
                print(snp_name, snps_to_include)
                assert False
            idx = varnames.index(snp_name)
            stored_indexes.append(idx)
            if not hardcalls:
                dosages = snp_dosages_or_hardcalls[:, 1] + 2*snp_dosages_or_hardcalls[:, 2]
                gt_temp[:, (snp_count - 1) % n_temp] = dosages
            else:
                gt_temp[:, (snp_count - 1) % n_temp] = snp_dosages_or_hardcalls
            assert not var_inclusion[snp_name]
            var_inclusion[snp_name] = True
            if snp_count % n_temp == 0:
                assert len(stored_indexes) == n_temp
                store_values(gts_dset, stored_indexes, gt_temp, covars, pheno_fname)
                gt_temp = np.full((n_samples, n_temp), np.nan)
                stored_indexes = []

        if len(stored_indexes) > 0:
            assert np.all(np.isnan(gt_temp[:, len(stored_indexes):]))
            gt_temp = gt_temp[:, :len(stored_indexes)]
            store_values(gts_dset, stored_indexes, gt_temp, covars, pheno_fname)

        print(f"Time: {time.time() - start}s")
        if sum(var_inclusion.values()) != len(var_inclusion):
            print(var_inclusion)
            print(snps_to_include)
            assert False


 
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('readme')
    parser.add_argument('gts_file')
    parser.add_argument('str_vcf')
    parser.add_argument('snp_bgen')
    parser.add_argument('varname_file', help='varnames should be first entry of each row when split on spaces')
    parser.add_argument('all_samples_file')
    parser.add_argument('phenotype_samples_file')
    parser.add_argument('phenotype')
    parser.add_argument('chrom')
    parser.add_argument('start')
    parser.add_argument('end')
    parser.add_argument('--varname-header', action='store_true', default=False, help='if true, skip the first line of the varname file')
    parser.add_argument('--pheno-fname')
    parser.add_argument('--shared-covar-fname')
    parser.add_argument('--pheno-out', help='where to save the pheno residuals')
    parser.add_argument('--hardcalls', action='store_true')
    args = parser.parse_args()

    assert bool(args.pheno_fname) == bool(args.shared_covar_fname) == bool(args.pheno_out)

    with tempfile.TemporaryDirectory(prefix='finemapping_load_gts', dir='.') as tempdir:
        print(f'tempdir name: {tempdir}', flush=True)
        readme = tempdir + '/' + args.readme.split('/')[-1]
        gts_file = tempdir + '/' + args.gts_file.split('/')[-1]
        if args.pheno_out:
            pheno_out = tempdir + '/' + args.pheno_out.split('/')[-1]
        else:
            pheno_out = None
        load_gts(
            readme,
            gts_file,
            args.str_vcf,
            args.snp_bgen,
            args.varname_file,
            args.all_samples_file,
            args.phenotype_samples_file,
            args.phenotype,
            args.chrom,
            args.start,
            args.end,
            args.varname_header,
            args.pheno_fname,
            args.shared_covar_fname,
            pheno_out,
            args.hardcalls
        )
        shutil.move(readme, args.readme)
        shutil.move(gts_file, args.gts_file)
        if args.pheno_out:
            shutil.move(pheno_out, args.pheno_out)

