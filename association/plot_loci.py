import argparse
import os
import sys

import cyvcf2
import matplotlib.pyplot as plt
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

    df = pd.read_csv(
        f'{assoc_dir}/results/{args.phenotype}.txt',
        usecols=['id', f'{args.phenotype}_residual'],
        sep=' '
    )

    with open(f'{assoc_dir}/phenotypes.txt') as units_file:
        first = True
        for line in units_file:
            if first:
                first = False
                continue
            this_phen, this_unit = line.split(":")
            if this_phen == args.phenotype:
                unit = this_unit

    fig, axs = plt.subplots((np.ceiling(len(args.locus)/2), 2),
                            figsize=(10, 5*(len(args.locus)/2))
                           )
    if len(axs) == 1:
        axs = np.array([[axs]])
    fig.set_title(
        f"Linear associations at individual loci against trait {args.phenotype} "
        "(covariate associations not plotted)"
    )

    for plot_idx, locus in enumerate(args.locus):
        ax = axs[plot_idx % 2, plot_idx/2]
        chrom, pos = locus.split(":")

        with open(f'{assoc_dir}/results/{args.phenotype}.txt') as results:
            coeff = None
            for line in results:
                line_chrom, line_pos = line.split()[:2]
                if chrom == line_chrom and pos == line_pos:
                    coeff = line.split()[5]
                    coeff = float(coeff)
                    break
            if coeff is None:
                raise ValueError(f"Couldn't find {locus} in results file"
                                 f" {results.name}!")
                

        vcf = cyvcf2.VCF(
            f'{ukb}/str_imputed/runs/{args.imputation_run_name}/vcfs/annotated_strs/chr{chrom}.vcf.gz'
        )
        record = next(iter(vcf(locus)))
        ref = record.REF
        period = record.INFO['PERIOD']
        period = int(period)
        bpdiffs = record.INFO['BPDIFFS']
        bpdiffs_repeat_units = bpdiffs/period

        idx_gts = record.genotype.array()[:, :-1]
        len_gts = np.zeros(idx_gts.shape)
        for idx, diff in enumerate(bpdiffs_repeat_units):
            len_gts[idx_gts == idx] = diff
        avg_len_gts = np.sum(len_gts, axis=1)
        '''
        command = f"""
        source ~/.bashrc ; 
        conda activate ukb ;
        bcftools query -f "REF %BPDIFFS %PERIOD" \\
            {ukb}/snpstr/info_field/chr{chrom}.vcf.gz
        """
        out = sp.run(command, shell=True, capture_output=True, text=True)
        if out.stderr:
            print("Command failed with stderr:\n", out.stderr)
            sys.exit(1)
        ref, bpdiffs, period = out.stdout.split()
        bpdiffs = bpdiffs.split(',')
        '''

        ax.set_title(f"Locus {locus}")
        ax.set_ylabel(unit)
        ax.set_xlabel("Diff from ref in repeat units")
        ax.annotate(
            f"Repeat period: {period}\tRef allele: {ref}\tLen: {len(ref)}",
            (0, -20)
        )

        phenotypes_by_genotypes = []
        axis_vals = []
        unique_gts = np.unique(avg_len_gts)
        for gt in unique_gts:
            axis_vals.append(gt)
            samples = vcf.samples[avg_len_gts == gt]
            phenotypes_by_genotypes.append(
                df.loc[samples, f'{args.phenotype}_residual']
            )

        ax.violin_plot(
            phenotypes_by_genotypes,
            axis_vals
        )
        ax.plot(unique_gts, unique_gts*coeff,
                label=f"Best fit line (slope == {coeff:.2e} {unit}/repeat unit")
        ax.legend()
    outloc = f'{assoc_dir}/loci_plots_{"_".join(args.locus)}.png'
    outloc = "_".join(outloc.split(':'))
    plt.savefig(outloc)


if __name__ == "__main__":
    main()
