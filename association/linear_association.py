#!/bin/env python3

"""
Run linear association of phenotypes
against covariates and STR genotypes

Reads the files using cyvcf2
and pumps them in to scipy
"""

import argparse
import datetime
import os
import os.path
import sys
import time

import cyvcf2
import dask
import dask.distributed
import dask_jobqueue
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
    df = pd.concat((ids_and_sex, pcs, const_df), axis=1)
    print("done")
    return df, indep_col_names


def load_height(df, readme, phenotypes):
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
    phenotypes.write("height:cm\n")
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
    # taller than tallest person or shorter than shortest adult
    filter_extremes = np.logical_or(height_df['height'] > 274,
                                    height_df['height'] < 54)
    n_filtered = np.sum(filter_extremes)
    readme.write(f"Filtering {n_filtered} height values that are taller than "
                 "the world's tallest person or shorter than the shortest.\n")
    height_df.loc[filter_extremes, 'height'] = np.nan

    try:
        return pd.merge(
            df,
            height_df,
            how="left",
            on="id"
        )
    finally:
        print("done")


def load_bilirubin(df, readme, phenotypes):
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
    phenotypes.write("total_bilirubin:umol/L\n")
    phenotypes.write("direct_bilirubin:umol/L\n")
    phenotypes.write("indirect_bilirubin:umol/L\n")
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

    filter_neg = np.sum(
        bilirubin_df.loc[:, ['total_bilirubin', 'direct_bilirubin', 'indirect_bilirubin']] < 0,
        axis=1
    )
    readme.write(f'Filtering {np.sum(filter_neg)} bilirubin values where '
                 'total or direct bilirubin was negative, or direct bilirubin '
                 'was greater than total bilirubin.')
    bilirubin_df.loc[filter_neg,
           ['total_bilirubin', 'direct_bilirubin', 'indirect_bilirubin']] = np.nan
    # max bilirubin value right now is 144. I don't know enough to say that
    # this is too high, so no max filter

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


@dask.delayed
def perform_association_subset(assoc_dir,
                               imputation_run_name,
                               region,
                               dep_var,
                               df):
    # make sure to use a local copy of df
    # so that dask doesn't try to change the copy we were passed
    df = df.copy(deep=True)
    chrom, poses = region.split(':')
    start, end = poses.split('-')
    region_string = f'{chrom}_{start}_{end}'
    with open(f'{assoc_dir}/run_logs/{dep_var}_{region_string}.log', 'w') as log, \
            open(f'{assoc_dir}/results/batches/{dep_var}_{region_string}.txt', 'w') as results:
        results.write("chrom pos #samples_GP_filtered #filtered_rare_alleles"
                      f" p_{dep_var} coeff_{dep_var} coeff_intercept\n")
        results.flush()
        vcf_floc = (f'{ukb}/str_imputed/runs/{imputation_run_name}/'
                    f'vcfs/strs_only/chr{chrom}.vcf.gz')
        vcf = cyvcf2.VCF(vcf_floc)

        vcf_region = vcf(region)
        n_loci = 0
        batch_time = 0
        batch_size = 50
        total_time = 0
        for locus in vcf_region:
            # setup for this locus
            n_loci += 1
            start_time = time.time()
            results.write(f"{locus.CHROM} {locus.POS}")

            # gt entries are sequence allele indexes
            idx_gts = locus.genotype.array()[:, :2]
            seq_allele_lens = [len(allele) for allele in locus.ALT]
            seq_allele_lens.insert(0, len(locus.REF))

            # get Beagle genotype probability
            aps = []
            for copy in 1, 2:
                ap = locus.format(f'AP{copy}')
                ref_prob = np.maximum(0, 1 - np.sum(ap, axis=1))
                ref_prob = ref_prob.reshape(-1, 1)
                ap = np.concatenate((ref_prob, ap), axis=1)
                ap = ap[np.arange(ap.shape[0]), idx_gts[:, (copy - 1)]]
                aps.append(ap)
            gp = np.multiply(aps[0], aps[1])

            # filter calls whose phased genotype probability
            # as estimated by Beagle is below .9
            # convert to float so we can use np.nan
            len_gts = np.zeros(idx_gts.shape, dtype=float)
            filtered_samples = gp < .9
            len_gts[filtered_samples, :] = np.nan
            n_filtered_samples = np.sum(filtered_samples)
            results.write(f" {n_filtered_samples}")

            # modify gt entries to be length alleles
            for seq_allele_idx, seq_allele_len in enumerate(seq_allele_lens):
                len_gts[idx_gts == seq_allele_idx] = seq_allele_len

            # filter length alleles with too few occurrences
            len_alleles, counts = np.unique(
                len_gts[~np.isnan(len_gts)],
                return_counts=True
            )
            filtered_rare_alleles = 0
            for len_allele_idx, len_allele in enumerate(len_alleles):
                if (counts[len_allele_idx] > 0 and
                        counts[len_allele_idx] < 5):
                    len_gts[len_gts == len_allele] = np.nan
                    filtered_rare_alleles += 1
            results.write(f" {filtered_rare_alleles}")

            # finalzie gts for regression
            avg_len_gt = np.sum(len_gts, axis=1)/2
            df['gt'] = avg_len_gt

            # make sure this locus doesn't have only one nonfiltered genotype
            not_nans = avg_len_gt[~np.isnan(avg_len_gt)]
            trivial_data = np.all(not_nans == not_nans[0])
            if trivial_data:
                # OLS thinks models with trivial data are infinitely
                # significant. For our purposes, they are not
                # significant at all.
                results.write(f" 1 nan")
                continue

            #do da regression
            model = OLS(
                df[f'{dep_var}_residual'],
                df[['gt', 'intercept']],
                missing='drop'
            )
            reg_result = model.fit()
            pval = reg_result.pvalues['gt']
            coef = reg_result.params['gt']
            intercept_coef = reg_result.params['intercept']
            results.write(f" {pval:.2e} {coef} {intercept_coef}")

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
    return region_string


@dask.delayed
def concat_batches(results_dir, dep_var, region_strings):
    with open(f'{results_dir}/{dep_var}.txt', 'w') as outfile:
        first = True
        for string in region_strings:
            with open(f'{results_dir}/batches/{dep_var}_{string}.txt') as infile:
                first_line = True
                for line in infile:
                    if first_line:
                        first_line = False
                        if first:
                            first = False
                            outfile.write(line)
                        continue
                    outfile.write(line)


def run_associations(assoc_dir, imputation_run_name, df, dep_vars):
    dask_output_dir = f'{assoc_dir}/dask_output_logs'
    os.mkdir(dask_output_dir)
    os.mkdir(f'{assoc_dir}/run_logs')
    os.mkdir(f'{assoc_dir}/results')
    os.mkdir(f'{assoc_dir}/results/batches')
    cluster = dask_jobqueue.PBSCluster(
        queue="condo",
        name="linear_associations",
        walltime="4:00:00",
        log_directory=dask_output_dir
    )
    cluster.adapt(maximum_jobs=300)
    client = dask.distributed.Client(cluster)

    print("Writing out results for each phenotype, region pair. "
          "See run_logs/<phenotype>_<region>.log for logs and "
          "See results/batches/<phenotype>_<region>.txt for results.",
         flush=True)

    # prep the data frame for being distributed across the cluster
    # see
    # https://docs.dask.org/en/latest/delayed-best-practices.html#avoid-repeatedly-putting-large-inputs-into-delayed-calls
    df = dask.delayed(df)

    #get chrom lengths
    chr_lens = np.loadtxt(
        f'{ukb}/misc_data/genome/chr_lens.txt',
        usecols=[1],
        skiprows=1
    )
    overall_dep_tasks = []
    for dep_var in dep_vars:
        tasks = []
        for chrom in (2, 3):
            for start_pos in range(1, int(chr_lens[chrom-1]), int(1e7)):
                #10Mbp regions
                tasks.append(perform_association_subset(
                    assoc_dir,
                    imputation_run_name,
                    f'{chrom}:{start_pos}-{start_pos + int(1e7) - 1}',
                    dep_var,
                    df
                ))
        overall_dep_tasks.append(concat_batches(
            f'{assoc_dir}/results', dep_var, tasks
        ))
    futures = client.compute(overall_dep_tasks)
    for future in futures:
        future.result() #block till all code is done executing
        # see here
        # https://distributed.dask.org/en/latest/manage-computation.html


def main():  # noqa: D103
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_filtering_run_name')
    parser.add_argument('imputation_run_name')
    parser.add_argument(
        'association_run_name',
        help="The name to give this association run"
    )
    """
    parser.add_argument(
        'regions',
        help=('comma separated list with elements of the form chr:start-end '
              'where chrom is mandatory but the rest are not')
    )
    """
    args = parser.parse_args()

    sample_filtering_run_name = args.sample_filtering_run_name
    imputation_run_name = args.imputation_run_name
    association_run_name = args.association_run_name
    # regions = args.regions.split(',')

    assoc_dir = f'{ukb}/association/runs/{association_run_name}'
    if os.path.exists(assoc_dir):
        print(f"Association with run name {association_run_name} already "
              f"exists!", file=sys.stderr)
        sys.exit(1)

    os.mkdir(assoc_dir)

    with open(f"{assoc_dir}/README", 'w') as readme, \
            open(f"{assoc_dir}/phenotypes.txt", 'w') as phenotypes:

        readme.write(f"Run args: {args}\n")
        today = datetime.datetime.now().strftime("%Y_%m_%d")
        readme.write(f"Run date: {today}\n")
        phenotypes.write("phenotype_name:unit\n")

        df, indep_vars = load_covars(readme)
        df = load_height(df, readme, phenotypes)
        df = load_bilirubin(df, readme, phenotypes)

        # Already calculated the largest unrelated subset of a filtering list
        # based on height. Will need to do this more intelligently when we
        # don't use the same samples for bilirubin and height
        floc = f'{ukb}/sample_qc/runs/{sample_filtering_run_name}/combined_unrelated.sample'
        readme.write(f"Subsetting to samples from file {floc}.\n")
        sample_subset = pd.read_csv(
            floc,
            usecols=[0],
            names=['id'],
            dtype=int,
            skiprows=1,
            sep=" "
        )
        df = pd.merge(df, sample_subset, how='inner', on='id')

        # order the observations in df
        # based on the order of samples in chr1 vcf
        print("Setting up samples ... ", end="", flush=True)
        vcf = cyvcf2.VCF(
            f'{ukb}/str_imputed/runs/{imputation_run_name}/'
            f'vcfs/strs_only/chr1.vcf.gz'
        )
        vcf_samples = list(sample.split("_")[0] for sample in
                           vcf.samples)
        samples_df = pd.DataFrame(
            vcf_samples,
            columns=['id'],
            dtype='int'
        )
        df = pd.merge(samples_df, df, how='left', on='id')
        df = df.astype({"id": "category"})
        print("done")

        # Calculate the number of samples used for each phenotype
        for phen in ('height', 'total_bilirubin'):
            readme.write(f'Number of samples for phenotype {phen}: ')
            readme.write(
                str(np.sum(np.all(
                    ~np.isnan(df.loc[:, [*indep_vars, phen]]),
                    axis=1)))
            )
            readme.write('\n')

        readme.write("Performing linear association against listed covariates"
                     " to get phenotypes residuals\n")
        # covariate columns:
        dep_vars = [
            'height',
            'total_bilirubin',
            #'direct_bilirubin',
            #'indirect_bilirubin'
        ]

        print("Running associations of phenotypes against covariates to get"
              " residuals ... ", end='', flush=True)
        for dep_var in dep_vars:
            with open(f'{assoc_dir}/{dep_var}_residual_results.txt', 'w') as res_results:
                model = OLS(
                    df[dep_var],
                    df[indep_vars],
                    missing='drop'
                )
                reg = model.fit()
                res_results.write("indep_var p coeff\n")
                for var in indep_vars:
                    res_results.write(f'{var} {reg.pvalues[var]:.2e} '
                                      f'{reg.params[var]}\n')
                predictions = reg.predict(df[indep_vars])
                df[f'{dep_var}_residual'] = df[dep_var] - predictions
        #Write out the csv we're going to be running associations on
        df.to_csv(
            f'{assoc_dir}/covars_and_phenotypes.txt',
            sep=' ',
            index=False
        )
        print("done", flush=True)

        readme.write(f"Using str calls from imputed run {imputation_run_name}\n")
        readme.write(
            "Filtering calls whose phased genotype probability as estimated "
            "by Beagle is less than 0.9. Number filtered this way is the "
            "#samples_GP_filtered column in filtering/<region>.txt"
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
        readme.write("Perofrming linear association of residual phenotypes"
                     "against allele length averaged over both haploid calls\n")
        readme.flush()

        run_associations(assoc_dir, imputation_run_name, df, dep_vars)


if __name__ == "__main__":
    main()

