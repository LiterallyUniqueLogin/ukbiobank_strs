#!/usr/bin/env python3

import argparse

import cyvcf2

parser = argparse.ArgumentParser()
parser.add_argument('vcf')
args = parser.parse_args()

vcf = cyvcf2.VCF(args.vcf)
n_samples = len(vcf.samples)

for record in vcf:
    allele_lens = [len(allele) for allele in [record.REF] + record.ALT]
    probs = {len_ : np.zeros(n_samples) for len_ in allele_lens}
