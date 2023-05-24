#!/usr/bin/env python3

import argparse

import cyvcf2
import numpy as np

import sample_utils

parser = argparse.ArgumentParser()
parser.add_argument('vcf')
parser.add_argument('sample_file')
parser.add_argument('all_samples_file')

args = parser.parse_args()

sample_idx = sample_utils.get_samples_idx(args.all_samples_file, args.sample_file)

vcf = cyvcf2.VCF(args.vcf)

print('chrom\tpos\tref_len\talleles\tcounts\tn')

for var in vcf:

    len_alleles = [len(allele) for allele in ([var.REF] + var.ALT)]

    allele_total_dosages = {
        _len: 0 for _len in np.unique(len_alleles)
    }

    for p in (1, 2):
        ap = var.format(f'AP{p}')
        allele_total_dosages[len_alleles[0]] += \
                np.sum(np.maximum(0, 1 - np.sum(ap[sample_idx, :], axis=1)))
        for i in range(ap.shape[1]):
            allele_total_dosages[len_alleles[i+1]] += np.sum(ap[sample_idx, i])

    print(
        f'{var.CHROM}\t'
        f'{var.POS}\t' +
        f'{len_alleles[0]}\t' +
        ','.join(str(len_) for len_ in np.unique(len_alleles)) + '\t' +
        ','.join(str(allele_total_dosages[allele]) for allele in np.unique(len_alleles)) + '\t' +
        str(np.sum(sample_idx))
    )

