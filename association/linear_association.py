#!/bin/env python3

"""
Run linear association of phenotypes
against covariates and STR genotypes

Reads the files using cyvcf2
and pumps them in to scipy
"""

import argparse
import datetime
import functools
import os
import os.path
import sys
import time
from typing import List, Optional

import cyvcf2
import dask
import dask.distributed
import dask_jobqueue
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import numpy.random
import sklearn.model_selection
import sklearn.neighbors
from statsmodels.regression.linear_model import OLS

import load_and_filter_genotypes

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

    # Age will be included as a covariate for all dependent variables
    # to prevent confounding.
    # If not, then if age was correlated with the dependent variable,
    # any loci which caused people to live longer or shorter would be
    # spuriously associated with the dependent variable.
    # However, each dependent variable measured during an assessment may
    # have been measured at one of multiple assessments.
    # As such, each measured dependent variable needs to specify which assessments
    # the age should be drawn from for each participant.
    # So we load ages for each assessment here but do not mark them as
    # independent variables to be directly included in the analysis.
    age_file_name = f'{ukb}/main_dataset/extracted_data/assessment_age.csv'
    readme.write(
        f"Adding ages at assessments. File: {age_file_name}\n"
    )
    assessment_age = pd.read_csv(
        age_file_name,
        names=["id", "age_assess_init", "age_assess_repeat", "age_assess_image"],
        dtype={"id": int,
               "age_assess_init": float,
               "age_assess_repeat": float,
               "age_assess_iamge": float},
        skiprows=1
    )
    df = pd.merge(df, assessment_age, how="inner", on="id")

    print("done")
    return df, indep_col_names


def load_height(df, readme, phenotypes):
    """
    Append height columns to the df.

    Height is measured in cm.

    Returns
    -------
    pd.DataFrame :
       Same as the input, but with three rows appended.
       The first is 'height' as a float
       (heights are measured in half cms)
       The second is 'height_sampling'
       and is categorical 0, 1 or 2
       specifying whether height was measured on the
       initial assessment, the first repeat assessment or
       the imaging assessment, respectively.
       Only first visit where height was retrieved is
       recorded.
       The last is 'height_age', namely the age of the ppt
       at the assessment specified by 'height_sampling'
    """
    print("Loading height ... ", end="", flush=True)
    phenotypes.write("height:cm\n")
    floc = f'{ukb}/main_dataset/extracted_data/height.txt'
    readme.write(
        f"Adding height phenotype, sampling visit and age at sampling visit. File: {floc}. "
        f"Only reporting first height measurement taken even if there "
        f"were multiple at different vists\n"
    )
    height_df = pd.read_csv(
        floc,
        names=["id", "height", "height_sampling"],
        dtype={"id": "int",
               "height": float,
               "height_sampling": float},
        skiprows=4,
        sep=" "
    )
    height_df.loc[pd.isna(height_df['height']), 'height_sampling'] = np.nan
    height_df.loc[pd.isna(height_df['height_sampling']), 'height'] = np.nan

    # taller than tallest person or shorter than shortest adult
    filter_extremes = np.logical_or(height_df['height'] > 274,
                                    height_df['height'] < 54)
    n_filtered = np.sum(filter_extremes)
    readme.write(f"Filtering {n_filtered} height values that are taller than "
                 "the world's tallest person or shorter than the shortest.\n")
    height_df.loc[filter_extremes, ['height', 'height_sampling']] = np.nan

    df = pd.merge(
        df,
        height_df,
        how="left",
        on="id"
    )

    readme.write("Adding age at height measurement as a covariate\n")
    start_idx = df.columns.get_loc('age_assess_init')
    has_height = ~np.isnan(df['height'])
    df.loc[has_height, 'height_age'] = \
        df.values[has_height.values, start_idx + df.loc[has_height, 'height_sampling'].astype(int)]
    print("done")
    return df, ['height_age']


def load_bilirubin(df, readme, phenotypes):
    """
    Append bilirubin columns to the df.

    Measured in umol/L

    Returns
    -------
    pd.DataFrame :
       Same as input df, but with four rows appended.
       The first recorded measurement is 'total_bilirubin'.
       The next column is 'total_bilirubin_age' and is the age
       at which that measurement was taken. The next column is
       'direct_bilirubin' which is set only for those ppts which
       had that measurement made at the same assessment as the total_bilirubin
       measurement. For those ppts, indirect bilirubin is calculated as
       total_bilirubin - direct_bilirubin.

       NOTE: it seems (table 1 :
           http://biobank.ctsu.ox.ac.uk/crystal/crystal/docs/serum_biochemistry.pdf
       ) that direct bilirubin values < 1 were filtered due to lack of
       precision, meaning that extant direct bilirubin levels are highly
       subject to bias.
    """
    print("Loading bilirubin ... ", end="", flush=True)
    phenotypes.write("total_bilirubin:umol/L\n")
    phenotypes.write("direct_bilirubin:umol/L\n")
    phenotypes.write("indirect_bilirubin:umol/L\n")
    floc = f'{ukb}/main_dataset/extracted_data/bilirubin.csv'
    readme.write(
        f"Adding direct, total and indirect_bilirubin phenotypes "
        f"and date of measurement, File: {floc}. "
        f"total_bilirubin is taken from the first assessment where it was "
        f"sampled, if any. direct_bilirubin is taken from that same assessment "
        f"if taken there. direct will not be taken from the other assessment. "
        f"direct will not be recorded if total is not."
        f"There are 70k total_bilirubin measurements with no corresponding "
        f"direct_bilirubin measurement and I do not know why, concerning."
        f"Almost all direct_bilirubin measurements have corresponding "
        f"total_bilirubin measurements, so no issue here."
        f"The indirect phenotype is simply calculated as total - direct.\n"
    )
    bilirubin_df = pd.read_csv(
        floc,
        names=["id", "dbil0", "dbil1", "tbil0", "tbil1"],
        header=0,
        dtype={"id": int,
               "dbil0" : float,
               "dbil1" : float,
               "tbil0" : float,
               "tbil1" : float},
        quotechar='"'
    )
    has_tbil0 = ~np.isnan(bilirubin_df.loc[:, 'tbil0'])
    bilirubin_df.loc[has_tbil0, 'bilirubin_assessment'] = 0
    bilirubin_df.loc[~has_tbil0, 'bilirubin_assessment'] = 1
    bilirubin_df.loc[has_tbil0, 'total_bilirubin'] = bilirubin_df.loc[has_tbil0, 'tbil0']
    bilirubin_df.loc[has_tbil0, 'direct_bilirubin'] = bilirubin_df.loc[has_tbil0, 'dbil0']
    bilirubin_df.loc[~has_tbil0, 'total_bilirubin'] = bilirubin_df.loc[~has_tbil0, 'tbil1']
    bilirubin_df.loc[~has_tbil0, 'direct_bilirubin'] = bilirubin_df.loc[~has_tbil0, 'dbil1']
    has_tbil = ~np.isnan(bilirubin_df.loc[:, 'total_bilirubin'])
    bilirubin_df.loc[~has_tbil, ['direct_bilirubin', 'bilirubin_assessment']] = np.nan
    bilirubin_df.drop(['tbil0', 'tbil1', 'dbil0', 'dbil1'], axis=1)

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

    df = pd.merge(
        df,
        bilirubin_df,
        how="left",
        on="id"
    )
    
    readme.write("Adding age at bilirubin measurement as a covariate\n")
    df.loc[has_tbil0, 'bilirubin_age'] = df.loc[has_tbil0, 'age_assess_init']
    df.loc[~has_tbil0, 'bilirubin_age'] = df.loc[~has_tbil0, 'age_assess_repeat']
    df.loc[~has_tbil, 'bilirubin_age'] = np.nan

    print("done")
    return df, ['bilirubin_age']


# copied from plot_stats branch of trtools
def PlotKDE(data: np.ndarray,
            xlabel: str,
            title: str,
            fname: str,
            strata_labels: Optional[List[str]] = None,
            random_state: int = 13):
    """
    Plots a kernel density estimation of the distribution.
    This is a smoother representation of the distribution
    than a historgram. Kernel bandwidth (which determins plot
    smoothness) is determined by cross validation (if more
    than 1000 loci, training is done on a subset of 1000
    chosen at random).
    Parameters
    ----------
    data:
        Either a 1D array of statistics to create a histogram of,
        or a 2D array where each column represents a different stratification
        of the data that will get a different line in the plot
    xlabel:
        the x label for the graph
    title:
        the title for the graph
    fname:
        the file name to save the graph. Must include the extension,
        and one that matplotlib will recognize so that it produces
        a file of that type.
    strata_labels:
        if data is 2D, then a different label for each column in data
    randome_state:
        used to control the splitting in the cross validation.
    """
    if len(data.shape) == 1:
        data = data.reshape(-1, 1)
    elif strata_labels is not None:
        assert len(strata_labels) == data.shape[1]

    n_strata = data.shape[1]

    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Probability density")

    # Fit and plot each stratum individually
    for col in range(n_strata):
        stratum = data[:, col]
        stratum = stratum[~np.isnan(stratum)]
        # only train on up to 1k loci for speed. Select a random subset
        # it would be nice to train on the same subset of loci for each strata
        # but that isn't feasible if some strata have many more nan's than
        # others
        max_loci = int(1e3)
        if len(stratum) > max_loci:
            rng = numpy.random.default_rng(random_state)
            stratum = stratum[rng.choice(len(stratum), size=max_loci, replace=False)]

        # code from
        # https://jakevdp.github.io/PythonDataScienceHandbook/05.13-kernel-density-estimation.html
        stratum = stratum.reshape(-1, 1)

        # Use gridsearch to choose the bandwidth
        kfold = sklearn.model_selection.KFold(
            n_splits=5,
            random_state=random_state,
            shuffle=True
        )
        bandwidths = 10 ** np.linspace(-1.8, 1.8, 20)
        grid = sklearn.model_selection.GridSearchCV(
            sklearn.neighbors.KernelDensity(kernel='gaussian'),
            {'bandwidth': bandwidths},
            cv=kfold
        )
        grid.fit(stratum)
        bandwidth = grid.best_params_['bandwidth']

        #compute the kde with the best bandwidth
        stratum = data[:, col] # now use all the data
        stratum = stratum[~np.isnan(stratum)]
        stratum = stratum.reshape(-1, 1)
        kde = sklearn.neighbors.KernelDensity(kernel='gaussian',
                                              bandwidth=bandwidth)
        kde.fit(stratum)
        min_val = np.min(stratum)
        max_val = np.max(stratum)
        eps = (max_val - min_val)/10e3
        xs = np.arange(min_val - eps, max_val + eps, eps)
        curve = np.exp(kde.score_samples(xs.reshape(-1, 1)))

        # plot
        if n_strata == 1:
            ax.fill_between(xs, curve)
        else:
            ax.plot(xs, curve, label=strata_labels[col])
    if n_strata > 1:
        ax.legend()
    plt.savefig(fname)


def plot_phenotype_by_sex(df, phenotype, fname, unit, title):
    print(f"Plotting {title} ... ", flush=True, end="")
    plot_df = df.loc[:, [phenotype, 'sex']].copy()
    plot_df.loc[:, f'male_{phenotype}'] = plot_df[phenotype]
    plot_df.loc[:, f'female_{phenotype}'] = plot_df[phenotype]
    plot_df.loc[plot_df.loc[:, 'sex'] == 2, f'male_{phenotype}'] = np.nan
    plot_df.loc[plot_df.loc[:, 'sex'] == 1, f'female_{phenotype}'] = np.nan
    PlotKDE(
        plot_df.loc[:, [f'male_{phenotype}', f'female_{phenotype}']].values,
        f'{phenotype} ({unit})',
        title,
        fname,
        ['male', 'female']
    )
    print("done", flush=True)


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

        n_loci = 0
        batch_time = 0
        batch_size = 50
        total_time = 0
        start_time = time.time()
        for len_gts, chrom, pos, n_filtered_samples, filtered_rare_alleles \
                in load_and_filter_genotypes.filtered_strs(
                    imputation_run_name, region
                ):
            n_loci += 1
            results.write(f"{chrom} {pos} {n_filtered_samples} "
                          f"{filtered_rare_alleles}")

            avg_len_gt = np.sum(len_gts, axis=1)/2
            df['gt'] = avg_len_gt

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
            start_time = time.time()
    return region_string


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

def run_associations_few_loci(assoc_dir, imputation_run_name, df, dep_vars,
                              loci):
    os.mkdir(f'{assoc_dir}/run_logs')
    os.mkdir(f'{assoc_dir}/results')
    os.mkdir(f'{assoc_dir}/results/batches')
   
    print("Writing out results for each phenotype, region pair. "
          "see run_logs/<phenotype>_<region>.log for logs and "
          "see results/batches/<phenotype>_<region>.txt for results.",
          flush=True)
    for dep_var in dep_vars:
        locus_strings = []
        for locus in loci:
            perform_association_subset(
                assoc_dir,
                imputation_run_name,
                locus,
                dep_var,
                df
            )
            locus_strings.append(locus.replace(':', '_').replace('-', '_'))
        concat_batches(f'{assoc_dir}/results', dep_var, locus_strings)


# see here:
# https://stackoverflow.com/questions/57118226/how-to-properly-use-dasks-upload-file-to-pass-local-code-to-workers
def _worker_upload(dask_worker, *, data, fname):
    dask_worker.loop.add_callback(
        callback=dask_worker.upload_file,
        comm=None,  # not used
        filename=fname,
        data=data,
        load=True)


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
    # currently dask cannot handle python files that are imported
    # locally. Have to upload it to all workers on the cluster,
    # but also have to make sure new workers
    # get the file too
    # see here:
    # https://stackoverflow.com/questions/57118226/how-to-properly-use-dasks-upload-file-to-pass-local-code-to-workers

    fname = 'load_and_filter_genotypes.py'
    with open(fname, 'rb') as f:
        data = f.read()

    client.register_worker_callbacks(
        setup=functools.partial(
            _worker_upload, data=data, fname=fname,
        )
    )


    print("Writing out results for each phenotype, region pair. "
          "see run_logs/<phenotype>_<region>.log for logs and "
          "see results/batches/<phenotype>_<region>.txt for results.",
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
        for chrom in range(1, 23):
            for start_pos in range(1, int(chr_lens[chrom-1]), int(1e7)):
                #10Mbp regions
                tasks.append(dask.delayed(perform_association_subset)(
                    assoc_dir,
                    imputation_run_name,
                    f'{chrom}:{start_pos}-{start_pos + int(1e7) - 1}',
                    dep_var,
                    df
                ))
        overall_dep_tasks.append(dask.delayed(concat_batches)(
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
    parser.add_argument(
        '--loci',
        help='comma separated list with elements of the form chr:pos'
    )
    args = parser.parse_args()

    sample_filtering_run_name = args.sample_filtering_run_name
    imputation_run_name = args.imputation_run_name
    association_run_name = args.association_run_name
    if args.loci:
        pre_format_loci = args.loci.split(',')
        loci = []
        for locus in pre_format_loci:
            chrom, pos = locus.split(':')
            loci.append(f'{chrom}:{pos}-{pos}')
    else:
        loci = None

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
        phenotypes.write("phenotype:unit\n")

        df, indep_vars = load_covars(readme)
        df, height_indep_vars = load_height(df, readme, phenotypes)
        df, bilirubin_indep_vars = load_bilirubin(df, readme, phenotypes)

        # covariate columns:
        dep_vars = [
            'height',
            'total_bilirubin',
            #'direct_bilirubin',
            #'indirect_bilirubin'
        ]
        pheno_indep_vars = {
            'height': height_indep_vars,
            'total_bilirubin': bilirubin_indep_vars
        }

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

        plot_dirs = {}
        for dep_var in dep_vars:
            plot_dirs[dep_var] = f'{assoc_dir}/{dep_var}_summary_plots'
            os.mkdir(plot_dirs[dep_var])

        plot_phenotype_by_sex(
            df,
            'total_bilirubin',
            f'{plot_dirs["total_bilirubin"]}/measured.png',
            'umol/L',
            'Measured total_bilirubin'
        )

        readme.write("Log transforming total_bilirubin\n")
        df.loc[:, 'total_bilirubin'] = np.log(df.loc[:, 'total_bilirubin'])

        plot_phenotype_by_sex(
            df,
            'total_bilirubin',
            f'{plot_dirs["total_bilirubin"]}/log_transform.png',
            'log(umol/L)',
            'Log transform total_bilirubin'
        )

        plot_phenotype_by_sex(
            df,
            'height',
            f'{plot_dirs["height"]}/measured.png',
            'cm',
            'Measured height'
        )

        readme.write("Performing linear association against listed covariates"
                     " to get phenotypes residuals\n")

        print("Running associations of phenotypes against covariates to get"
              " residuals ... ", end='', flush=True)
        units = {'height': 'cm', 'total_bilirubin': 'log(umol/L)'}
        for dep_var in dep_vars:
            curr_indep_vars = indep_vars.copy()
            curr_indep_vars.extend(pheno_indep_vars[dep_var])
            with open(f'{assoc_dir}/{dep_var}_residual_results.txt', 'w') as res_results:
                model = OLS(
                    df[dep_var],
                    df[curr_indep_vars],
                    missing='drop'
                )
                reg = model.fit()
                res_results.write("indep_var p coeff\n")
                for var in curr_indep_vars:
                    res_results.write(f'{var} {reg.pvalues[var]:.2e} '
                                      f'{reg.params[var]}\n')
                predictions = reg.predict(df[curr_indep_vars])
                df[f'{dep_var}_residual'] = df[dep_var] - predictions
        #Write out the csv we're going to be running associations on
        df.to_csv(
            f'{assoc_dir}/covars_and_phenotypes.txt',
            sep=' ',
            index=False
        )
        print("done", flush=True)
        for dep_var in dep_vars:
            plot_phenotype_by_sex(
                df,
                f'{dep_var}_residual',
                f'{plot_dirs[dep_var]}/residual.png',
                f'residual {units[dep_var]}',
                f'Residual {dep_var}'
            )

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
        readme.write("Performing linear association of residual phenotypes "
                     "against allele length averaged over both haploid calls."
                     " Coefficients are reported as change in phenotype per "
                     "difference from reference allele measured in repeat "
                     "copies.\n")
        readme.flush()

        if loci is not None:
            run_associations_few_loci(assoc_dir, imputation_run_name, df, dep_vars, loci)
        else:
            run_associations(assoc_dir, imputation_run_name, df, dep_vars)

if __name__ == "__main__":
    main()

