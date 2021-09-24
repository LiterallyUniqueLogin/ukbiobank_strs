import argparse
import os

import numpy as np

import load_and_filter_genotypes as lfg
import python_array_utils as utils

ukb = os.environ['UKB']

def get_str_dosages(str_dosages_dict):
    return np.sum([_len*np.sum(dosages, axis=1) for
                   _len, dosages in str_dosages_dict.items()], axis=0)

def load_gts(imputation_run_name, phenotype, chrom, start, end):
    with open(f'{ukb}/sample_qc/runs/{phenotype}/combined_unrelated.sample') as samples_file:
        samples = np.array([line.strip() for line in samples_file][1:], dtype=int).reshape(-1, 1)

    imp_snp_samples_filepath = f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample'
    with open(imp_snp_samples_filepath) as imp_snp_samples_file:
        imp_snp_samples = np.array([line.split()[0] for line in imp_snp_samples_file][2:], dtype=int)

    samples_indicator = np.concatenate((samples, samples), axis=1)
    samples_merge = utils.merge_arrays(imp_snp_samples, samples_indicator)
    assert samples_merge.shape[1] == 2
    sample_idx = ~np.isnan(samples_merge[:, 1])

    var_names = []
    gts = [] # also includes other covars

    pheno_covars = np.load(f'{ukb}/traits/subset_transformed_phenotypes/{phenotype}.npy')
    pheno_covar_names = [line.strip() for line in open(f'{ukb}/traits/phenotypes/{phenotype}_covar_names.txt')]

    # reorder samples to the proper order
    pheno_covars = utils.merge_arrays(samples_merge[sample_idx, 0:1], pheno_covars)
    
    var_names.extend(pheno_covar_names)
    gts.extend(pheno_covars[:, i] for i in range(2, pheno_covars.shape[1]))

    region = f'{chrom}:{start}-{end}'
    str_itr = lfg.load_strs(imputation_run_name, region, sample_idx)
    next(str_itr) # skip the locus details list

    for str_dosages_dict, _, _, str_pos, str_locus_filtered, _ in str_itr:
        if str_locus_filtered:
            continue
        # TODO skip if p-value is too high?
        var_names.append(f'STR_{str_pos}')
        gts.append(get_str_dosages(str_dosages_dict))

    snp_itr = lfg.load_imputed_snps(region, sample_idx)

    for snp_dosages, alleles, _, snp_pos, snp_filtered, _ in snp_itr:
        if snp_filtered:
            continue
        var_names.append(f'SNP_{snp_pos}_' + '_'.join(alleles))
        # TODO skip if p-value is too high?
        gts.append(snp_dosages[:, 1] + 2*snp_dosages[:, 2])

    outdir = f'{ukb}/finemapping/susie_results/{phenotype}/{chrom}_{start}_{end}'
    outfname = f'{outdir}/gts.tab'
    np.savefile(outfname, np.stack(gts), delimiter='\t')
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
