#!/usr/bin/env python3

import os

import argparse
import cyvcf2
import numpy as np

ukb = os.environ['UKB']

# total len, len polyA, len TA, len CA, len T(A/G), SNP category
allele_details = np.array([
[52, 2, 18, 14, 18, 0],
[48, 2, 14, 16, 16, 0],
[48, 2, 14, 14, 18, 0],
[49, 5, 12, 14, 18, 0],
[49, 3, 14, 16, 16, 0],
[49, 3, 14, 14, 18, 0],
[49, 1, 16, 14, 18, 0],
[50, 6, 12, 16, 16, 0],
[50, 6, 12, 14, 18, 0],
[50, 4, 14, 16, 16, 0],
[50, 4, 14, 14, 18, 0],
[50, 4, 14, 14, 18, 1],
[50, 2, 16, 16, 16, 0],
[50, 2, 16, 14, 18, 0],
[51, 5, 14, 14, 18, 0],
[51, 3, 16, 14, 18, 2],
[51, 3, 16, 14, 18, 0],
[51, 1, 16, 16, 18, 0],
[51, 1, 18, 14, 18, 0],
[52, 2, 18, 14, 18, 0],
[53, 5, 14, 16, 18, 0],
[53, 5, 16, 14, 18, 0],
[53, 3, 16, 16, 18, 0],
[53, 3, 18, 14, 18, 3],
[53, 3, 18, 14, 18, 4],
[53, 3, 18, 14, 18, 5],
[53, 3, 18, 12, 20, 0],
[53, 1, 18, 16, 18, 0],
[53, 1, 20, 14, 18, 0],
[54, 6, 16, 14, 18, 0],
[54, 4, 18, 14, 18, 0],
[54, 2, 20, 14, 18, 0],
[55, 5, 18, 14, 18, 0],
[55, 3, 18, 16, 18, 0],
[55, 3, 20, 14, 18, 0],
])
varnames = [
    'PACSIN2_complex_repeat_len',
    'PACSIN2_polyA_STR_partial_len',
    'PACSIN2_TA_STR_len',
    'PACSIN2_CA_STR_len',
    'PACSIN2_T(A|G)_STR_partial_len',
    'PACSIN2_SNP_43385897_C_G',
    'PACSIN2_SNP_43385885_T_C',
    'PACSIN2_SNP_43385893_C_A_SNP_43385924_G_C',
    'PACSIN2_SNP_43385903_C_A',
    'PACSIN2_SNP_43385917_T_G'
]

assert np.all(np.sum(allele_details[:, 1:-1], axis=1) == allele_details[:, 0])

def load_dosages(var, gts):
    both_aps = []
    for chrom in range(1, 3):
        aps = var.format(f'AP{chrom}')
        ref_aps = np.maximum(0, 1 - np.sum(aps, axis=1))
        aps = np.hstack((ref_aps.reshape(-1, 1), aps))
        both_aps.append(aps)
        assert aps.shape[1] == allele_details.shape[0]
        for a in range(allele_details.shape[0]):
            for len_ in range(5):
                gts[:, len_] += allele_details[a, len_]*aps[:, a]
            snp = allele_details[a, -1]
            if snp != 0:
                gts[:, 4 + snp] += aps[:, a]
    both_aps = np.stack(both_aps, axis=1)
    np.save(f'{ukb}/finemapping/PACSIN2/idx_dosages.npy', both_aps)
    np.save(f'{ukb}/finemapping/PACSIN2/dosages.npy', gts)

def load_hardcalls(var, gts):
    alleles = var.genotype.array()[:, :-1]
    np.save(f'{ukb}/finemapping/PACSIN2/idx_hardcalls.npy', alleles)
    for chrom in range(2):
        for a in range(allele_details.shape[0]):
            samps = alleles[:, chrom] == a
            for len_ in range(5):
                gts[samps, len_] += allele_details[a, len_]
            snp = allele_details[a, -1]
            if snp != 0:
                gts[samps, 4 + snp] += 1
    np.save(f'{ukb}/finemapping/PACSIN2/hardcalls.npy', gts)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hardcalls', default=False, action='store_true')
    args = parser.parse_args()

    vcf = cyvcf2.VCF(f'{ukb}/str_imputed/runs/first_pass/vcfs/annotated_strs/chr22.vcf.gz')
    region = vcf('22:43385872')
    var = next(region)
    # first 5 are STR len gts, next 5 are SNP gts
    gts = np.zeros((len(vcf.samples), 10))
    if not args.hardcalls:
        load_dosages(var, gts)
    else:
        load_hardcalls(var, gts)
        
if __name__ == '__main__':
    main()
