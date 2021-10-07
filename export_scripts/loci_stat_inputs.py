#!/usr/bin/env python3

import argparse
import os
import sys

import cyvcf2
import numpy as np
import scipy.stats

ukb = os.environ['UKB']

sys.path.insert(0, f'{ukb}/../trtools/repo')

import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('str_imputation_run_name')
    parser.add_argument('chrom')
    args = parser.parse_args()
    str_imputation_run_name = args.str_imputation_run_name
    chrom = args.chrom

    with open(f'{ukb}/export_scripts/intermediate_results/chr{chrom}_loci_summary.tab', 'w') as output:
        main_helper(output, str_imputation_run_name, chrom)

def main_helper(output, str_imputation_run_name, chrom):
    vcf_fname = (f'{ukb}/str_imputed/runs/{str_imputation_run_name}/'
                             f'vcfs/annotated_strs/chr{chrom}.vcf.gz')
    vcf = cyvcf2.VCF(vcf_fname)

    subset_samples = []
    with open(f'{ukb}/sample_qc/runs/no_phenotype/combined.sample') as samp_file:
        next(samp_file)
        for line in samp_file:
            subset_samples.append(line.strip())

    all_samples = [l[0] for l in np.char.split(vcf.samples, '_')]
    samples = np.isin(all_samples, subset_samples)

    output.write('chr\tpos\tallele_dist\tentropy\theterozygosity\tmultiallelicness\n')
    for record in vcf:
        if record.INFO.get('PERIOD') is None:
            continue
        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype='beagle-hipstr')

        len_alleles = [trrecord.ref_allele_length] + trrecord.alt_allele_lengths
        total_subset_dosages = {
            len_: 0 for len_ in np.unique(len_alleles)
        }
        for p in (1, 2):
            ap = trrecord.format[f'AP{p}']
            total_subset_dosages[len_alleles[0]] += \
                    np.sum(np.maximum(0, 1 - np.sum(ap[samples, :], axis=1)))
            for i in range(ap.shape[1]):
                total_subset_dosages[len_alleles[i+1]] += np.sum(ap[samples, i])

        for len_ in total_subset_dosages:
            total_subset_dosages[len_] /= np.sum(samples)*2

        entropy = scipy.stats.entropy(list(total_subset_dosages.values()), base=2)
        heterozygosity = 1 - np.sum(val**2 for val in total_subset_dosages.values())
        multiallelicness = sum(sorted(total_subset_dosages.values())[:-2])

        output.write(
            trrecord.chrom + '\t' + str(trrecord.pos) + '\t' +
            str(total_subset_dosages) + '\t' + str(entropy) + '\t' +
            str(heterozygosity) + '\t' + str(multiallelicness) + '\n'
        )

if __name__ == '__main__':
    main()
