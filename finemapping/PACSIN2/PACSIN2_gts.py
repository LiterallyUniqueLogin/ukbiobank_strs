#!/usr/bin/env python3

import os

import argparse
import cyvcf2
import numpy as np

from PACSIN2_varnames import varnames

ukb = os.environ['UKB']

# total len, len polyA, len TA, len CA, len T(A/G), SNP category
allele_details = np.array([
[53, 3, 18, 14, 18, 0],
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


assert np.all(np.sum(allele_details[:, 1:-1], axis=1) == allele_details[:, 0])

def load_dosages(var, nsamples):
    both_aps = []
    varname_to_len_to_gts = {}
    for idx, varname in enumerate(varnames[:5]):
        varname_to_len_to_gts[varname] = {}
        for len_ in np.unique(allele_details[:, idx]):
            varname_to_len_to_gts[varname][len_] = np.zeros((nsamples, 2))
    for chrom in range(1, 3):
        aps = var.format(f'AP{chrom}')
        ref_aps = np.maximum(0, 1 - np.sum(aps, axis=1))
        aps = np.hstack((ref_aps.reshape(-1, 1), aps))
        both_aps.append(aps)
        assert aps.shape[1] == allele_details.shape[0]
        for idx, varname in enumerate(varnames[:5]):
            for len_ in varname_to_len_to_gts[varname]:
                varname_to_len_to_gts[varname][len_][:, chrom-1] += np.sum(aps[:, allele_details[:, idx] == len_], axis=1)
            '''
                gts[:, len_] += allele_details[a, len_]*aps[:, a]
            snp = allele_details[a, -1]
            if snp != 0:
                gts[:, 4 + snp] += aps[:, a]
            '''
    both_aps = np.stack(both_aps, axis=1)
    np.save(f'{ukb}/finemapping/PACSIN2/gts/idx_dosages.npy', both_aps)
    for varname in varname_to_len_to_gts:
        for len_ in varname_to_len_to_gts[varname]:
            np.save(f'{ukb}/finemapping/PACSIN2/gts/{varname}_{len_}.npy', varname_to_len_to_gts[varname][len_])

def load_hardcalls(var, nsamples):
    alleles = var.genotype.array()[:, :-1]
    np.save(f'{ukb}/finemapping/PACSIN2/gts/idx_hardcalls.npy', alleles)
    gts = np.zeros((nsamples, 10))
    for chrom in range(2):
        for a in range(allele_details.shape[0]):
            samps = alleles[:, chrom] == a
            for len_ in range(5):
                gts[samps, len_] += allele_details[a, len_]
            snp = allele_details[a, -1]
            if snp != 0:
                gts[samps, 4 + snp] += 1
    np.save(f'{ukb}/finemapping/PACSIN2/gts/hardcalls.npy', gts)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hardcalls', default=False, action='store_true')
    args = parser.parse_args()

    vcf = cyvcf2.VCF(f'{ukb}/str_imputed/runs/first_pass/vcfs/annotated_strs/chr22.vcf.gz')
    region = vcf('22:43385872')
    var = next(region)
    nsamples = len(vcf.samples)
    # first 5 are STR len gts, next 5 are SNP gts
    if not args.hardcalls:
        load_dosages(var, nsamples)
    else:
        load_hardcalls(var, nsamples)
        
if __name__ == '__main__':
    main()
