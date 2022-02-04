#!/usr/bin/env python3

import argparse
import os

import numpy as np

import load_and_filter_genotypes as lfg
import PACSIN2_varnames
import sample_utils
import test_PACSIN2

def process_STR(itr, pos):
    try:
        dosage_dict, unique_alleles, chrom, found_pos, locus_filtered, locus_details = next(itr)
    except StopIteration as si:
        raise ValueError(f"Did not find a STR at position {pos}") from si

    assert int(pos) == int(found_pos)
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

def PACSIN2_itr(pos):
    itr = test_PACSIN2.get_gt_itr(slice(None))
    next(itr)
    for dosage_dict, unique_alleles, chrom, found_pos, locus_filtered, locus_details in itr:
        if found_pos != pos:
            continue
        yield (dosage_dict, unique_alleles, chrom,  found_pos, locus_filtered, locus_details)

def main():
    ukb = os.environ['UKB']

    parser = argparse.ArgumentParser()
    parser.add_argument('outprefix')
    parser.add_argument('phenotype')
    parser.add_argument('chr')
    parser.add_argument('--imputation-run-name')
    parser.add_argument('--STRs', nargs='+', default=[])
    parser.add_argument('--imputed-SNPs', nargs='+', default=[])
    parser.add_argument('--PACSIN2-STRs', nargs='+', default=[])

    args = parser.parse_args()

    assert len(args.STRs) + len(args.imputed_SNPs) + len(args.PACSIN2_STRs) > 0

    if len(args.STRs) > 0:
        assert args.imputation_run_name

    if len(args.PACSIN2_STRs) > 0:
        assert args.chr == '22'

    samples_array = sample_utils.get_all_samples().reshape(-1).astype(float)
    assert len(samples_array) == 487409

    variant_names = []
    gts = []
    gts.append(samples_array)

    for STR in args.STRs:
        itr = lfg.load_strs(args.imputation_run_name, f'{args.chr}:{STR}-{STR}', slice(None))
        next(itr) # ignore names of locus_details
        this_gts, this_var_name = process_STR(itr, STR)
        variant_names.append(this_var_name)
        gts.append(this_gts)

    for STR in args.PACSIN2_STRs:
        itr = PACSIN2_itr(STR)
        this_gts, this_var_name = process_STR(itr, STR)
        variant_names.append(this_var_name)
        gts.append(this_gts)

    for iSNP in args.imputed_SNPs:
        pos, ref, alt = iSNP.split('_')
        itr = lfg.load_imputed_snps(f'{args.chr}:{pos}-{pos}', slice(None))
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
