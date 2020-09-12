import argparse
import os
import subprocess as sp
import sys

import cyvcf2
import numpy as np
import pandas as pd

ukb = os.environ['UKB']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("imputation_run_name")
    parser.add_argument("association_run_name")
    parser.add_argument("phenotype")
    parser.add_argument("locus", nargs="+",
                        help="should be represented as chrom:pos")

    args = parser.parse_args()

    assoc_dir = f'{ukb}/association/runs/{args.association_run_name}'
    plot_dir = f'{assoc_dir}/locus_plots/{args.phenotype}'
    res_col = f'{args.phenotype}_residual'

    os.makedirs(plot_dir, exist_ok=True)

    df = pd.read_csv(
        f'{assoc_dir}/covars_and_phenotypes.txt',
        usecols=['id', res_col],
        sep=' '
    )
    data = df.to_numpy()
    data = data[~np.isnan(data[:, 1]), :]

    with open(f'{assoc_dir}/phenotypes.txt') as units_file:
        first = True
        success = False
        for line in units_file:
            if first:
                first = False
                continue
            phen, unit = line.split(":")
            if phen == args.phenotype:
                with open(f'{plot_dir}/info.txt', 'w') as phenotype_info:
                    phenotype_info.write(f'{unit}\n')
                    success = True
                    break
        if not success:
            print(f"Couldn't find phenotype {args.phenotype}", file=sys.stderr)
            sys.exit(1)

    for plot_idx, locus in enumerate(args.locus):
        chrom, pos = locus.split(":")

        with open(f'{assoc_dir}/results/{args.phenotype}.txt') as results:
            coeff = None
            for line in results:
                line_chrom, line_pos = line.split()[:2]
                if chrom == line_chrom and pos == line_pos:
                    p, coeff, intercept = line.split()[4:7]
                    coeff = float(coeff)
                    intercept = float(intercept)
                    p = float(p)
                    break
            if coeff is None:
                raise ValueError(f"Couldn't find {locus} in results file"
                                 f" {results.name}!")

        vcf = cyvcf2.VCF(
            f'{ukb}/str_imputed/runs/{args.imputation_run_name}/vcfs/annotated_strs/chr{chrom}.vcf.gz'
        )
        record = next(iter(vcf(locus)))
        with open(f'{plot_dir}/{chrom}_{pos}.info', 'w') as plot_info:
            plot_info.write(f'{record.REF}\n')
            period = record.INFO["PERIOD"]
            plot_info.write(f'{period}\n')
            plot_info.write(f'{coeff}\n')
            plot_info.write(f'{intercept}\n')
            plot_info.write(f'{p}\n')

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

        gt_df = pd.DataFrame(data={
            'gt': avg_len_gts,
            'id': samples
        })

        merged_df = pd.merge(
            df,
            gt_df,
            how="inner",
            on="id"
        )

        merged_df.to_csv(
            f'{plot_dir}/{chrom}_{pos}.csv',
            columns=['gt', res_col],
            index=False
        )

        command = f"""
        source ~/.bashrc ;
        conda activate ukb ;
        Rscript {ukb}/association/plot_loci_helper.R \
                {plot_dir} {chrom} {pos} {args.phenotype}
        """
        sp.run(command, check=True, shell=True)

if __name__ == "__main__":
    main()
