#!/usr/bin/env python3

import argparse

import cyvcf2
import numpy as np

import trtools.utils.tr_harmonizer as trh 

parser = argparse.ArgumentParser()
parser.add_argument('vcf')
parser.add_argument('white_brits')
parser.add_argument('outdir')
args = parser.parse_args()

chrom = 'chr11'
pos = 119206290
vcf = cyvcf2.VCF(args.vcf)

samples = np.isin(vcf.samples, [line.split('_')[0] for line in open(args.white_brits).readlines()])
assert np.sum(samples) > 100 

var = next(vcf(f'{chrom}:{pos}'))
rec = trh.HarmonizeRecord('hipstr', var)
seqs = rec.GetStringGenotypes()[samples, :-1]

hom_pure_samples = np.array(vcf.samples)[samples][np.all(
    np.char.find(seqs, 'T', start=3) != 3,
    axis=1
)]
hom_impure_samples = np.array(vcf.samples)[samples][np.all(
    np.char.find(seqs, 'T', start=3) == 3,
    axis=1
)]

with open(f'{args.outdir}/cbl_pure_white_brits.samples', 'w') as pure_out:
    pure_out.write('ID\n')
    pure_out.write('\n'.join(hom_pure_samples))
    pure_out.write('\n')

with open(f'{args.outdir}/cbl_impure_white_brits.samples', 'w') as impure_out:
    impure_out.write('ID\n')
    impure_out.write('\n'.join(hom_impure_samples))
    impure_out.write('\n')

