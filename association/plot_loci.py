import argparse
import os
import subprocess as sp

import cyvcf2
import numpy as np
import pandas as pd

import python_array_utils as utils

ukb = os.environ['UKB']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("imputation_run_name")
    parser.add_argument("phenotype")
    parser.add_argument('transform_step', choices=['original', 'rin', 'adjusted'])
    parser.add_argument("locus", nargs="+",
                        help="should be represented as chrom:pos")

    args = parser.parse_args()
    phenotype = args.phenotype

    plot_dir = f'{ukb}/association/locus_plots/{phenotype}'
    os.makedirs(plot_dir, exist_ok=True)

    if args.transform_step == 'original':
        phenotypes = np.load(f'{ukb}/traits/phenotypes/{phenotype}.npy')[:, :2]
        with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt') as unit_file:
            unit = next(unit_file).strip()
    elif args.transform_step == 'rin':
        phenotypes = np.load(f'{ukb}/traits/subset_rin_phenotypes/{phenotype}.npy')[:, :2]
        unit = 'rank inverse normalized'
    elif args.transform_step == 'adjusted':
        phenotypes = np.load(f'{ukb}/traits/adjusted_srin_phenotypes/{phenotype}_linear.npy')[:, :2]
        unit = 'adjusted rank inverse normalized'

    for locus in args.locus:
        chrom, pos = locus.split(":")

        vcf = cyvcf2.VCF(
            f'{ukb}/str_imputed/runs/{args.imputation_run_name}/vcfs/annotated_strs/chr{chrom}.vcf.gz'
        )
        record = next(iter(vcf(locus)))
        period = record.INFO["PERIOD"]

        bpdiffs = list(record.INFO['BPDIFFS'])
        bpdiffs.insert(0, 0)
        bpdiffs = np.array(bpdiffs)
        bpdiffs_repeat_units = bpdiffs/period

        idx_gts = record.genotype.array()[:, :-1]
        len_gts = np.zeros(idx_gts.shape)
        for idx, diff in enumerate(bpdiffs_repeat_units):
            len_gts[idx_gts == idx] = diff
        avg_len_gts = np.sum(len_gts, axis=1)/2

        samples = np.array(vcf.samples)
        samples = np.char.partition(samples, '_')[:, 0]
        samples = samples.astype(int)

        gts_array = np.concatenate((samples.reshape(-1, 1), avg_len_gts.reshape(-1, 1)), axis=1)

        outarray = utils.merge_arrays(phenotypes, gts_array)

        np.savetxt(
            f'{plot_dir}/{chrom}_{pos}_{args.transform_step}.csv',
            outarray[:, 1:],
            delimiter=','
        )

        command = f"""
        source ~/.bashrc ;
        conda activate ukb ;
        Rscript {ukb}/association/plot_loci_helper.R \
                {plot_dir} {chrom} {pos} {args.phenotype} {args.transform_step} {unit} {period} {record.REF}
        """
        sp.run(command, check=True, shell=True)

if __name__ == "__main__":
    main()
