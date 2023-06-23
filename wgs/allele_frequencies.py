#!/usr/bin/env python3

import argparse
import json

import cyvcf2
import numpy as np
import trtools.utils.tr_harmonizer as trh

parser = argparse.ArgumentParser()
parser.add_argument('vcf')
parser.add_argument('sample_lists', nargs=6)
args = parser.parse_args()

itrs = []
for sample_list in args.sample_lists:
    samples = [line.split()[0] for line in open(sample_list).readlines() if line.strip()]
    vcf = cyvcf2.VCF(args.vcf, samples=samples)
    itrs.append(trh.TRRecordHarmonizer(vcf))

print('\t'.join(
    ['chrom', 'pos_hg38',
     'white_brit_n_called_samples', 'white_brit_freqs',
     'black_n_called_samples,' 'black_freqs',
     'south_asian_n_called_samples', 'south_asian_freqs',
     'chinese_n_called_samples', 'chinese_freqs',
     'irish_n_called_samples', 'irish_freqs',
     'white_other_n_called_samples', 'white_other_freqs']
))

done = None
n_recs = 0
while done is not True:
    done = None
    for itr_count, itr in enumerate(itrs):
        try:
            rec = next(itr)
        except StopIteration:
            assert done is not False
            done = True
            continue
        done = False
        if itr_count == 0:
            print(f'{rec.chrom}\t{rec.pos}', end='')
        print(f'\t{np.sum(rec.GetCalledSamples())}\t{json.dumps(rec.GetAlleleFreqs())}', end='')
    print()
    n_recs += 1
    if n_recs == 10:
        break

