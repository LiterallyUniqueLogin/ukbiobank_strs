#!/usr/bin/env python3

import argparse
import os
import time

import numpy as np
from statsmodels.stats.weightstats import DescrStatsW

import load_and_filter_genotypes as lfg
import python_array_utils as utils

ukb = os.environ['UKB']

parser = argparse.ArgumentParser()
parser.add_argument('chrom')
args = parser.parse_args()
chrom = args.chrom

with open(f'{ukb}/sample_qc/runs/white_brits/no_phenotype/combined_unrelated.sample') as samples_file:
    samples = np.array([line.strip() for line in samples_file][1:], dtype=int).reshape(-1, 1)

imp_snp_samples_filepath = f'{ukb}/array_imputed/ukb46122_imp_chr1_v3_s487283.sample'
with open(imp_snp_samples_filepath) as imp_snp_samples_file:
    imp_snp_samples = np.array([line.split()[0] for line in imp_snp_samples_file][2:], dtype=int).reshape(-1, 1)

samples_indicator = np.concatenate((samples, samples), axis=1)
samples_merge = utils.merge_arrays(imp_snp_samples, samples_indicator)
assert samples_merge.shape[1] == 2
sample_idx = ~np.isnan(samples_merge[:, 1])
n_samples = np.sum(sample_idx)

itr = lfg.load_strs('first_pass', f'{chrom}:1-1000000000', sample_idx, details=False)
outloc = f'{ukb}/side_analyses/length_confusion/chr{chrom}.tab'
start_time = time.time()
with open(f'{outloc}.temp', 'w') as out:
    out.write('pos\tchance_of_length_confusion\tavg_abs_length_confusion\tnormalized_avg_abs_length_confusion\n')
    for str_count, (str_dosages_dict, _, _, str_pos, str_locus_filtered, _) in enumerate(itr):
        if str_count % 100 == 0:
            print(f'Working on STR num {str_count+1} pos {str_pos} ', end ='')
            if str_count != 0:
                print(f'time/str: {(time.time()-start_time)/str_count:.2}s', flush=True)
            else:
                print('', flush=True)
            out.flush()
        if str_locus_filtered:
            continue
        chance_of_length_confusion = np.zeros((n_samples, 2))
        avg_length_confusion = np.zeros((n_samples, 2))
        lengths = np.zeros(n_samples*2*len(str_dosages_dict))
        all_probabilities = np.zeros(n_samples*2*len(str_dosages_dict))
        for allele_num, (len1, dosages1) in enumerate(str_dosages_dict.items()):
            lengths[(allele_num*n_samples*2):((allele_num+1)*n_samples*2)] = len1
            all_probabilities[(allele_num*n_samples*2):((allele_num+1)*n_samples*2)] = \
                dosages1.reshape(-1)
            for len2, dosages2 in str_dosages_dict.items():
                if len1 == len2:
                    continue
                chance_of_length_confusion += \
                    np.multiply(dosages1, dosages2)
                avg_length_confusion += \
                    np.multiply(dosages1, dosages2)*np.abs(len1-len2)
        std = DescrStatsW(lengths, weights=all_probabilities).std
        out.write(f'{str_pos}\t{np.mean(chance_of_length_confusion)}\t{np.mean(avg_length_confusion)}\t{np.mean(avg_length_confusion)/std}\n')

os.rename(f'{outloc}.temp', outloc)
