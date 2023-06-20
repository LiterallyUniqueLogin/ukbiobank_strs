#!/usr/bin/env python3

import argparse

import numpy as np

import load_and_filter_genotypes as lfg
import sample_utils

def process_STR(itr, pos):
    try:
        dosage_dict, unique_alleles, chrom, found_pos, locus_filtered, locus_details = next(itr)
    except StopIteration as si:
        raise ValueError(f"Did not find a STR at position {pos}") from si

    if locus_filtered:
        raise ValueError(f"STR at position {pos} was filtered for reason {locus_filtered}")

    try:
        next(itr)
    except StopIteration:
        pass
    else:
        raise ValueError(f"Expected only one STR at position {pos}")

    this_var_gts = np.sum(
        [_len*np.sum(dosages, axis=1) for _len, dosages in dosage_dict.items()],
        axis=0
    )
    this_var_gts = (this_var_gts - np.mean(this_var_gts))/np.std(this_var_gts)
    return (this_var_gts, f"STR_{pos}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outprefix')
    parser.add_argument('all_samples')
    parser.add_argument('chr')
    parser.add_argument('--str-vcf')
    parser.add_argument('--snp-bgen')
    parser.add_argument('--snp-mfi')
    parser.add_argument('--STRs', nargs='*', default=[]) # will be ignored if empty
    parser.add_argument('--imputed-SNPs', nargs='*', default=[]) # will be ignored if empty

    args = parser.parse_args()

    assert len(args.STRs) + len(args.imputed_SNPs)

    if len(args.STRs) > 0:
        assert args.str_vcf

    if len(args.imputed_SNPs) > 0:
        assert args.snp_bgen
        assert args.snp_mfi

    samples_array = sample_utils.get_samples(args.all_samples).reshape(-1).astype(float)
    assert len(samples_array) == 487409

    variant_names = []
    gts = []
    gts.append(samples_array)

    for STR in args.STRs:
        itr = lfg.load_strs(args.str_vcf, f'{args.chr}:{STR}-{STR}', slice(None))
        next(itr) # ignore names of locus_details
        this_gts, this_var_name = process_STR(itr, STR)
        variant_names.append(this_var_name)
        gts.append(this_gts)

    for iSNP in args.imputed_SNPs:
        pos, ref, alt = iSNP.split('_')
        itr = lfg.load_imputed_snps(args.snp_bgen, args.snp_mfi, f'{args.chr}:{pos}-{pos}', slice(None))
        next(itr) # ignore names of locus_details
        try:
            unique_alleles = []
            while not np.all(unique_alleles == [ref, alt]):
                dosages, unique_alleles, chrom, outpos, locus_filtered, locus_details = next(itr)
        except StopIteration as si:
            raise ValueError(f"Did not find imputed SNP {iSNP}") from si

        if locus_filtered:
            raise ValueError(f"imputed SNP {iSNP} was filtered for reason {locus_filtered}")
        assert outpos == int(pos)

        variant_names.append(f"imputed_SNP_{iSNP}")
        this_var_gts = dosages[:, 1] + 2*dosages[:, 2]
        this_var_gts = (this_var_gts - np.mean(this_var_gts))/np.std(this_var_gts)
        gts.append(this_var_gts)

    np.save(
        f'{args.outprefix}.npy',
        np.stack(gts, axis=1)
    )
    with open(f'{args.outprefix}_varnames.txt', 'w') as f:
        f.write('sample_ID ' + ' '.join(variant_names) + '\n')

    with open(f'{args.outprefix}_README.txt', 'w') as f:
        f.write('Computed single dosages per sample for STRs in number of repeats, single '
                'dosages per sample for SNPs as number of copies of the alternate allele.\n')

if __name__ == '__main__':
    main()
