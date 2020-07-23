#!/bin/env python3

"""
A simple for loop.

Reads the files using cyvcf2
and pumps them in to scipy
"""

import argparse
import os
import os.path
import sys

import cyvcf2
import pandas as pd
import numpy as np

ukb = os.environ['UKB']


def load_covars(readme):
    """
    Load the sex and population PC covariates.
    TODO include age

    Returns
    -------
    pd.DataFrame :
        A categorical column 'id' which serves as the logical index
        An integer column 'sex' with M = 1 and F = 2
        Float columns pc1 ... pc40 for the population structure pcs

    Notes
    -----
    If interested in the source of this data, read the READMEs in the
    directories of the loaded files.
    """
    print("Loading covariates ... ", end="", flush=True)
    floc = f'{ukb}/microarray/ukb46122_cal_chr1_v2_s488282.fam'
    cols = (0, 4)
    readme.write(
        f"Adding participant ID index and sex covariate to df. File: {floc}, cols: {cols}\n"
    )
    ids_and_sex = pd.read_csv(
        floc,
        usecols=cols,
        header=None,
        names=['id', 'sex'],
        dtype=int,
        sep=" "
    )

    floc = f'{ukb}/misc_data/EGA/ukb_sqc_v2.txt'
    cols = list(range(25, 65))
    readme.write(
        f"Adding PC covariates 1-40. File: {floc}, cols: {cols}. Participants in same "
        f"order as previous file.\n"
    )
    pcs = pd.read_csv(
        floc,
        sep=" ",
        header=None,
        names=list(f"pc{col}" for col in range(1, 41)),
        usecols=cols,
        dtype=float
    ) #156MB

    try:
        return pd.concat((ids_and_sex, pcs), axis=1)
    finally:
        print("done")


def load_height(df, readme):
    """
    Append height columns to the df.

    Height is measured in cm.

    Returns
    -------
    pd.DataFrame :
       Same as the input, but with two rows appended.
       The first is 'height' as a float
       (heights are measured in half cms)
       The second is 'height_sampling'
       and is categorical 0, 1 or 2
       specifying which visit height was retrieved on.
       Only first visit where height was retrieved is
       recorded.
    """
    print("Loading height ... ", end="", flush=True)
    floc = f'{ukb}/main_dataset/extracted_data/height.txt'
    readme.write(
        f"Adding height phenotype and sampling visit. File: {floc}. "
        f"Only reporting first height measurement taken even if there "
        f"were multiple at different vists\n"
    )
    height_df = pd.read_csv(
        floc,
        names=["id", "height", "height_sampling"],
        dtype={"id": "int",
               "height": float,
               "height_sampling": "category"},
        skiprows=4,
        sep=" "
    )

    try:
        return pd.merge(
            df,
            height_df,
            how="left",
            on="id"
        )
    finally:
        print("done")


def load_bilirubin(df, readme):
    """
    Append bilirubin columns to the df.

    Measured in umol/L

    Returns
    -------
    pd.DataFrame :
       Same as input df, but with two float rows appended.
       They are 'direct bilirubin' and
       'total bilirubin'. 'indirect bilirubin'
       can be calculated as total - direct.
       These only reflect measurements made at the first
       visit, follow up measurements are ignored
       as I do not know if bilirubin values change as people
       age or if I should be adjusting for that.
       Unsure why someone would get a value for
       total or direct bilirubin but not the other,
       so for now make sure both exist or NaN them both
    """
    print("Loading bilirubin ... ", end="", flush=True)
    floc = f'{ukb}/main_dataset/extracted_data/bilirubin.csv'
    readme.write(
        f"Adding direct and total bilirubin phenotypes. File: {floc}. "
        f"Only recording first visit's value, regardless of whether "
        f"it is not present and/or there is a follow up visit value. "
        f"NaNing total values where direct is missing (or visa versa) "
        f"as an overly conservative approach to making sure we have clean "
        f"data as a I don't know why we'd be missing one but not both. "
        f"This effectively cuts down 70k total bilirubin measurements "
        f"which don't have direct counterparts\n"
    )
    bilirubin_df = pd.read_csv(
        floc,
        names=["id", "direct bilirubin", "total bilirubin"],
        header=0,
        usecols=[0, 1, 3],
        dtype={"id": "int",
               "direct bilirubin" : float,
               "total bilirubin" : float},
        quotechar='"'
    )
    direct_nan = np.isnan(bilirubin_df['direct bilirubin'])
    total_nan = np.isnan(bilirubin_df['total bilirubin'])
    bilirubin_df.loc[direct_nan, 'total bilirubin'] = np.nan
    bilirubin_df.loc[total_nan, 'direct bilirubin'] = np.nan

    """
    both_nan = np.logical_and(direct_nan, total_nan)
    print("Direct nan count", np.sum(direct_nan))  # 103892
    print("total nan count", np.sum(total_nan))  # 34945
    print("both nan count", np.sum(both_nan))  # 34601
    """

    try:
        return pd.merge(
            df,
            bilirubin_df,
            how="left",
            on="id"
        )
    finally:
        print("done")


def main():  # noqa: D103
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_filtering_run_name')
    parser.add_argument('imputation_run_name')
    parser.add_argument(
        'association_run_name',
        help="The name to give this association run"
    )
    args = parser.parse_args()

    sample_filtering_run_name = args.sample_filtering_run_name
    imputation_run_name = args.imputation_run_name
    association_run_name = args.association_run_name

    assoc_dir = f'{ukb}/association/runs/{association_run_name}'
    if os.path.exists(assoc_dir):
        print(f"Association with run name {association_run_name} already "
               "exists!", file=sys.stderr)
        sys.exit(1)

    os.mkdir(assoc_dir)

    with open(f"{assoc_dir}/REAMDE", 'w') as readme, \
            open(f"{assoc_dir}/run.log", 'w') as log, \
            open(f"{assoc_dir}/results.txt", 'w') as results:
        df = load_covars(readme)
        df = load_height(df, readme)
        df = load_bilirubin(df, readme)

        """
        height_nan = np.isnan(df['height'])
        bilirubin_nan = np.isnan(df['direct bilirubin'])
        both_nan = np.logical_and(height_nan, bilirubin_nan)
        print("height nan count", np.sum(height_nan))  # 1487
        print("bilirubin nan count", np.sum(bilirubin_nan))  # 93578
        print("both nan count", np.sum(both_nan))  # 437
        """

        # For ease of use, only work with samples which are
        # present for both
        # Fix later
        readme.write(
            "Removing samples which are missing height or bilirubin values "
            "so we can run both together. "
            "This cuts down on ~100k ppts who have height but not bilirubin. "
            "Will do this more intelligently later.\n"
        )
        df = df[~np.isnan(df['height'])]
        df = df[~np.isnan(df['total bilirubin'])]

        # Already calculated the largest unrelated subset of a filtering list
        # based on height. Will need to do this more intelligently when we
        # don't use the same samples for bilirubin and height
        floc = f'{ukb}/sample_qc/runs/{sample_filtering_run_name}/combined_unrelated.sample'
        readme.write(f"Subsetting to samples from file {floc}. Left with "
                     "~275k samples\n")
        sample_subset = pd.read_csv(
            floc,
            usecols=[0],
            names=['id'],
            dtype=int,
            skiprows=1,
            sep=" "
        )
        df = pd.merge(df, sample_subset, how='inner', on='id')

        df = df.astype({"id": "category"})

        # Run the association assuming all phenotypes are present simultaneously
        # QQ plot and Manhattan plot
        # How can a permutation test (for a QQ plot) still be above the line?

        readme.write(f"Using str calls from imputed run {imputation_run_name}\n")
        readme.wrtie(
            "Filtering calls whose phased genotype probability as estimated "
            "by Beagle is less than 0.9. Number filtered this way is the "
            "#samples_GP_filtered column in results.txt"
        )
        readme.write(
            "Convert and collapse sequence alleles to length alleles. "
            "Removing length alleles with 4 or fewer occurrences. "
            "For a single association, removing samples with any removed "
            "alleles. Then convert sample genotypes to their average "
            "allele length. #filtered_rare_alleles in results.txt is the "
            "number of length alleles with between 1 and 4 copies in all "
            "samples which were filtered due to rarity.\n"
        )
        readme.write("Performing linear association against listed covariates"
                     " and average allele length\n")

        #Write results header
        results.write("chrom pos #samples_GP_filtered #filtered_rare_alleles")
        # 2 for bilirubin, 3 for height
        #TODO for sex chroms take into account ploidy
        for chrom in 2, 3:
            vcf_floc = (f'{ukb}/str_imputed/runs/{imputation_run_name}/'
                        f'vcfs/strs_only/chr{chrom}.vcf.gz')
            for locus in cyvcf2.VCF(vcf_floc):
                results.write(f"{locus.CHROM} {locus.POS}")
                # gt entries are sequence allele indexes
                # convert to float so we can use np.nan
                gt = locus.genotype.array()[:2].astype(float)
                seq_allele_lens = [len(allele) for allele in locus.ALT]
                seq_allele_lens.insert(0, len(locus.REF))

                # filter calls whose phased genotype probability
                # as estimated by Beagle is below .9
                aps = []
                for copy in 1, 2:
                    ap = locus.format(f'AP{copy}')
                    ref_prob = np.maximum(0, 1 - np.sum(ap, axis=1))
                    ap = np.concatenate((ref_prob, ap), axis=1)
                    ap = ap[:, gt[:, (copy - 1)]]
                    aps.append(ap)
                gp = np.mult(aps[0], aps[1])

                filtered_samples = gp < .9
                gt[filtered_samples, :] = np.nan
                n_filtered_samples = np.sum(filtered_samples)
                results.write(f" {n_filtered_samples}")

                # modify gt entries to be length alleles
                for seq_allele_idx, seq_allele_len in enumerate(seq_allele_lens):
                    gt[gt == seq_allele_idx] = seq_allele_len

                # filter length alleles with too few occurrences
                len_alleles, counts = np.unique(gt, return_counts=True)[1]
                filtered_rare_alleles = 0
                for len_allele_idx, len_allele in enumerate(len_alleles):
                    if (counts[len_allele_idx] > 0 and
                            counts[len_allele_idx] < 5):
                        gt[gt == len_allele] = np.nan
                        filtered_rare_alleles += 1
                results.write(f" {filtered_rare_alleles}")

                avg_len = np.sum(gt, axis=1)/2
                results.write("\n")


if __name__ == "__main__":
    main()

