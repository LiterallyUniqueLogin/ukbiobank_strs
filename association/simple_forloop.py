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
import time

import cyvcf2
import pandas as pd
import numpy as np
from statsmodels.regression.linear_model import OLS

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
        and an 'intercept' column of ones

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
    pc_col_names = list(f"pc{col}" for col in range(1, 41))
    readme.write(
        f"Adding PC covariates 1-40. File: {floc}, cols: {cols}. Participants in same "
        f"order as previous file.\n"
    )
    pcs = pd.read_csv(
        floc,
        sep=" ",
        header=None,
        names=pc_col_names,
        usecols=cols,
        dtype=float
    ) #156MB

    const = np.ones((len(pcs.index), 1))
    const_df = pd.DataFrame(const, columns=["intercept"], dtype="int16")
    indep_col_names = ['intercept', 'sex']
    indep_col_names.extend(pc_col_names)
    try:
        return (
            pd.concat((ids_and_sex, pcs, const_df), axis=1),
            indep_col_names
        )
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
       They are 'direct_bilirubin' and
       'total_bilirubin'. 'indirect_bilirubin'
       can be calculated as total - direct.
       These only reflect measurements made at the first
       visit, follow up measurements are ignored
       as I do not know if bilirubin values change as people
       age or if I should be adjusting for that.
       Unsure why someone would get a value for
       total or direct_bilirubin but not the other,
       so for now make sure both exist or NaN them both
    """
    print("Loading bilirubin ... ", end="", flush=True)
    floc = f'{ukb}/main_dataset/extracted_data/bilirubin.csv'
    readme.write(
        f"Adding direct, total and indirect_bilirubin phenotypes. File: {floc}. "
        f"Only recording first visit's value, regardless of whether "
        f"it is not present and/or there is a follow up visit value. "
        f"NaNing total values where direct is missing (or visa versa) "
        f"as an overly conservative approach to making sure we have clean "
        f"data as a I don't know why we'd be missing one but not both. "
        f"This effectively cuts down 70k total_bilirubin measurements "
        f"which don't have direct counterparts. The indirect phenotype"
        f" is simply calculated as total - direct.\n"
    )
    bilirubin_df = pd.read_csv(
        floc,
        names=["id", "direct_bilirubin", "total_bilirubin"],
        header=0,
        usecols=[0, 1, 3],
        dtype={"id": "int",
               "direct_bilirubin" : float,
               "total_bilirubin" : float},
        quotechar='"'
    )
    direct_nan = np.isnan(bilirubin_df['direct_bilirubin'])
    total_nan = np.isnan(bilirubin_df['total_bilirubin'])
    bilirubin_df.loc[direct_nan, 'total_bilirubin'] = np.nan
    bilirubin_df.loc[total_nan, 'direct_bilirubin'] = np.nan

    # Add indirect:
    bilirubin_df['indirect_bilirubin'] = (
        bilirubin_df['total_bilirubin'] -
        bilirubin_df['direct_bilirubin']
    )

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
    parser.add_argument(
        'regions',
        help=('comma separated list with elements of the form chr:start-end '
              'where chrom is mandatory but the rest are not')
    )
    args = parser.parse_args()

    sample_filtering_run_name = args.sample_filtering_run_name
    imputation_run_name = args.imputation_run_name
    association_run_name = args.association_run_name
    regions = args.regions.split(',')

    assoc_dir = f'{ukb}/association/runs/{association_run_name}'
    if os.path.exists(assoc_dir):
        print(f"Association with run name {association_run_name} already "
              f"exists!", file=sys.stderr)
        sys.exit(1)

    os.mkdir(assoc_dir)

    with open(f"{assoc_dir}/README", 'w') as readme, \
            open(f"{assoc_dir}/run.log", 'w') as log, \
            open(f"{assoc_dir}/results.txt", 'w') as results:

        readme.write(f"Run args: {args}")

        df, indep_vars = load_covars(readme)
        df = load_height(df, readme)
        df = load_bilirubin(df, readme)

        """
        height_nan = np.isnan(df['height'])
        bilirubin_nan = np.isnan(df['direct_bilirubin'])
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
        df = df[~np.isnan(df['total_bilirubin'])]

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
        readme.write(
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


        # covariate columns:
        dep_vars = [
            'height',
            'total_bilirubin',
            'direct_bilirubin',
            'indirect_bilirubin'
        ]

        #Write results header
        results.write("chrom pos #samples_GP_filtered #filtered_rare_alleles")
        for dep_var in dep_vars:
            results.write(f" p_{dep_var} coeff_{dep_var}")
        results.write("\n")
        results.flush()

        n_loci = 0
        batch_time = 0
        batch_size = 50
        total_time = 0
        # 2 for bilirubin, 3 for height
        #TODO for sex chroms take into account ploidy
        samples_set_up = False

        for region in regions:

            chrom, in_chrom_region = region.split(':')
            if '-' not in in_chrom_region:
                # cyvcf2 interprets 2:500
                # as chromosome 2, bp 500 and on
                # even though the way bcftools inteprets that is as
                # chrom 2, bp 500 exactly
                # so change the syntax we pass to cyvcf2 to get
                # it onboard
                region = f"{chrom}:{in_chrom_region}-{in_chrom_region}"

            vcf_floc = (f'{ukb}/str_imputed/runs/{imputation_run_name}/'
                        f'vcfs/strs_only/chr{chrom}.vcf.gz')
            vcf = cyvcf2.VCF(vcf_floc)

            # order the observations in df
            # based on the order of samples in the first vcf
            if not samples_set_up:
                print("Setting up samples ...", end="", flush=True)
                vcf_samples = list(sample.split("_")[0] for sample in
                                   vcf.samples)
                samples_df = pd.DataFrame(
                    vcf_samples,
                    columns=['id'],
                    dtype='int'
                )
                df = pd.merge(samples_df, df, how='left', on='id')

                indep_var_array = df[indep_vars].to_numpy()
                indep_vars.insert(0, 'gt')

                n_samples = np.sum(df['intercept'] == 1)
                readme.write(f"Running associations with {n_samples} "
                             f"samples\n")
                readme.flush()
                samples_set_up = True
                print("done")

            print(f"Writing out results for region {region}. See the log for "
                  f"updates.")

            vcf_region = vcf(region)
            for locus in vcf_region:
                # setup for this locus
                n_loci += 1
                start_time = time.time()
                results.write(f"{locus.CHROM} {locus.POS}")

                # gt entries are sequence allele indexes
                gt = locus.genotype.array()[:, :2]
                seq_allele_lens = [len(allele) for allele in locus.ALT]
                seq_allele_lens.insert(0, len(locus.REF))

                # get Beagle genotype probability
                aps = []
                for copy in 1, 2:
                    ap = locus.format(f'AP{copy}')
                    ref_prob = np.maximum(0, 1 - np.sum(ap, axis=1))
                    ref_prob = ref_prob.reshape(-1, 1)
                    ap = np.concatenate((ref_prob, ap), axis=1)
                    ap = ap[np.arange(ap.shape[0]), gt[:, (copy - 1)]]
                    aps.append(ap)
                gp = np.multiply(aps[0], aps[1])

                # filter calls whose phased genotype probability
                # as estimated by Beagle is below .9
                # convert to float so we can use np.nan
                gt = gt.astype(float)
                filtered_samples = gp < .9
                gt[filtered_samples, :] = np.nan
                n_filtered_samples = np.sum(filtered_samples)
                results.write(f" {n_filtered_samples}")

                # modify gt entries to be length alleles
                for seq_allele_idx, seq_allele_len in enumerate(seq_allele_lens):
                    gt[gt == seq_allele_idx] = seq_allele_len

                # filter length alleles with too few occurrences
                len_alleles, counts = np.unique(
                    gt[~np.isnan(gt)],
                    return_counts=True
                )
                filtered_rare_alleles = 0
                for len_allele_idx, len_allele in enumerate(len_alleles):
                    if (counts[len_allele_idx] > 0 and
                            counts[len_allele_idx] < 5):
                        gt[gt == len_allele] = np.nan
                        filtered_rare_alleles += 1
                results.write(f" {filtered_rare_alleles}")

                # set gts for regression
                avg_len_gt = np.sum(gt, axis=1)/2

                #do da regression
                for dep_var in dep_vars:
                    model = OLS(
                        df[dep_var],
                        np.concatenate((avg_len_gt.reshape(-1, 1),
                                        indep_var_array),
                                       axis=1),
                        missing='drop'
                    )
                    model.exog_names[:] = indep_vars
                    reg_result = model.fit()
                    pval = reg_result.pvalues['gt']
                    coef = reg_result.params['gt']
                    results.write(f" {pval} {coef}")

                # output results
                results.write("\n")
                results.flush()
                duration = time.time() - start_time
                total_time += duration
                batch_time += duration
                if n_loci % batch_size == 0:
                    log.write(
                        f"time/locus (last {batch_size}): "
                        f"{batch_time/batch_size}s\n"
                        f"time/locus ({n_loci} total loci): {total_time/n_loci}s\n"
                    )
                    log.flush()
                    batch_time = 0


if __name__ == "__main__":
    main()

