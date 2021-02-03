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
import subprocess as sp
import sys
import time
import tracemalloc
from typing import List, Optional

import csaps
import cyvcf2
import dask
import dask.distributed
import dask_jobqueue
import matplotlib.pyplot as plt
import numpy as np
import numpy.random
import sklearn.model_selection
import sklearn.neighbors
from statsmodels.regression.linear_model import OLS

import load_and_filter_genotypes

ukb = os.environ['UKB']

def merge_arrays(a, b):
    # assume first column of each array is id, but not necessarily same order
    # goal is return do a left outer join b

    assert len(set(a[:, 0]).intersection(b[:,0])) > 1000
    assert len(set(a[:, 0])) == a.shape[0]
    assert len(set(b[:, 0])) == b.shape[0]

    b = b[np.isin(b[:, 0], a[:, 0])]
    matches = np.isin(a[:, 0], b[:, 0])

    a_sort = np.argsort(a[matches, 0])
    b_match_sorted = np.searchsorted(a[matches, 0], b[:, 0], sorter=a_sort)

    new_data = np.full((a.shape[0], b.shape[1] - 1), np.nan)
    new_data[matches, :] = b[np.argsort(b_match_sorted), 1:][np.argsort(a_sort), :]

    return np.concatenate((
        a,
        new_data
    ), axis=1)
    

def load_covars(readme):
    """
    Load the sex and population PC covariates.
    TODO include age

    Returns
    -------
    data : np.ndarray
        A 2D float array with one row per sample and columns in the following order:
        id, sex with M = 1 and F = 2, pc1 ... pc40 for the population
        structure pcs, and an 'intercept' column of ones
    colnames: List[str]
        The names of the columns in the returned array
    indep_cols : List[str]
        A boolean array denoting which of the returned columsn are independent
        variables (in this case, all but the id)

    Notes
    -----
    If interested in the source of this data, read the READMEs in the
    directories of the loaded files.
    """
    print("Loading covariates ... ", end="", flush=True)
    floc = f'{ukb}/microarray/ukb46122_cal_chr1_v2_s488282.fam'
    cols = (0, 4)
    readme.write(
        f"Loading participant ID index and sex covariate. File: {floc}, cols: {cols}\n"
    )
    readme.flush()
    ids_and_sex = np.genfromtxt(
        floc,
        usecols=cols,
        delimiter=" "
    )
    col_names = ['id', 'sex']

    floc = f'{ukb}/misc_data/EGA/ukb_sqc_v2.txt'
    cols = list(range(25, 65))
    col_names.extend(list(f"pc{col}" for col in range(1, 41)))
    readme.write(
        f"Adding PC covariates 1-40. File: {floc}, cols: {cols}. Participants in same "
        f"order as previous file.\n"
    )
    readme.flush()
    pcs = np.genfromtxt(
        floc,
        delimiter=" ",
        usecols=cols
    ) #156MB
    

    const_array = np.ones((pcs.shape[0], 1))
    col_names.append('intercept')

    # these arrays are all in the same row order, so just concatenate
    data = np.concatenate((ids_and_sex, pcs, const_array), axis=1)

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
    readme.flush()
    with open(age_file_name) as age_file:
        assessment_age = np.genfromtxt(
            (line.replace('"', '') for line in age_file),
            skip_header=1,
            delimiter=','
        )
    data = merge_arrays(data, assessment_age)
    col_names.extend(["age_assess_init", "age_assess_repeat",
                     "age_assess_image"])
    # don't include id or any of the ages
    indep_cols = list(range(1, len(col_names)-3))
    print("done")

    return data, col_names, indep_cols


def load_height(data, col_names, readme, phenotypes):
    """
    Append height columns to the array

    Height is measured in cm.

    Returns
    -------
    data : np.ndarray
       Same as the input, but with three columns appended.
       The first is 'height' (heights are measured in half cms)
       The second is 'height_sampling'
       and is categorical 0, 1 or 2
       specifying whether height was measured on the
       initial assessment, the first repeat assessment or
       the imaging assessment, respectively.
       Only first visit where height was retrieved is
       recorded.
       The last is 'height_age', namely the age of the ppt
       at the assessment specified by 'height_sampling'
    col_names : List[str]
        Same as input, but extended with 'height', 'height_sampling' and
        'height_age'
    height_indep_vars : List[str]
        ['height_age']
    """
    print("Loading height ... ", end="", flush=True)
    phenotypes.write("height:cm\n")
    floc = f'{ukb}/main_dataset/extracted_data/height.txt'
    readme.write(
        f"Adding height phenotype, sampling visit and age at sampling visit. File: {floc}. "
        f"Only reporting first height measurement taken even if there "
        f"were multiple at different vists\n"
    )
    readme.flush()
    # cols "id", "height", "height_sampling"
    height_data = np.genfromtxt(
        floc,
        skip_header=4,
        delimiter=" "
    )

    # taller than tallest person or shorter than shortest adult
    filter_extremes = np.logical_or(height_data[:, 1] > 274,
                                    height_data[:, 1] < 54)
    n_filtered = np.sum(filter_extremes)
    readme.write(f"Filtering {n_filtered} height values that are taller than "
                 "the world's tallest person or shorter than the shortest.\n")
    readme.flush()
    height_data[np.ix_(filter_extremes, [1,2])] = np.nan
    assert not np.any(np.isnan(height_data))

    data = merge_arrays(data, height_data)

    readme.write("Adding age at height measurement as a covariate\n")
    readme.flush()
    data = np.concatenate(
        (data, np.full((data.shape[0], 1), np.nan)),
        axis=1
    )
    col_names.extend(['height', 'height_sampling', 'height_age'])

    start_idx = col_names.index('age_assess_init')
    has_height = ~np.isnan(data[:, col_names.index('height')])
    data[has_height, -1] = data[
        has_height,
        start_idx + data[has_height, col_names.index('height_sampling')].astype(int)
    ]

    print("done")
    return data, col_names, ['height_age']


def load_bilirubin(data, col_names, readme, phenotypes):
    """
    Append bilirubin columns to the array

    Measured in umol/L

    Returns
    -------
    data : np.ndarray
       Same as input, but with four columns appended.
       The first two are tbil0 and tbil1 which are just intermediate
       calculations.
       The next is 'total_bilirubin'.
       The last is 'total_bilirubin_age' and is the age
       at which that measurement was taken.
    col_names : List[str]
        Same as input, but extended with ['tbil0', 'tbil1', 
        'total_bilirubin, 'total_bilirubin_age']
    bilirubin_indep_vars : List[str]
        ['total_bilirubin_age']
    """
    print("Loading bilirubin ... ", end="", flush=True)
    phenotypes.write("total_bilirubin:log(umol/L)\n")
    floc = f'{ukb}/main_dataset/extracted_data/bilirubin.csv'
    readme.write(
        f"Adding total_bilirubin phenotype "
        f"and date of measurement, File: {floc}. "
        f"total_bilirubin is taken from the first assessment where it was "
        f"sampled, if any.\n"
    )
    readme.flush()
    # cols id", "tbil0", "tbil1"
    with open(floc) as bilirubin_csv:
        bilirubin_data = np.genfromtxt(
            (line.replace('"', '') for line in bilirubin_csv),
            skip_header=1,
            usecols=[0,3,4],
            delimiter=","
        )

    data = merge_arrays(data, bilirubin_data)

    data = np.concatenate(
        (data, np.full((data.shape[0], 2), np.nan)),
        axis=1
    )
    col_names.extend(['tbil0', 'tbil1', 'total_bilirubin', 'bilirubin_age'])

    use_tbil0 = ~np.isnan(data[:, -4])
    use_tbil1 = ~use_tbil0 & ~np.isnan(data[:, -3])
    data[use_tbil0, -2] = data[use_tbil0, -4]
    data[use_tbil1, -2] = data[use_tbil1, -3]

    # max bilirubin value right now is 144. I don't know enough to say that
    # this is too high, so no max filter
    # I've checked, none are negative

    readme.write("Adding age at bilirubin measurement as a covariate\n")
    readme.flush()
    data[use_tbil0, -1] = data[use_tbil0, col_names.index('age_assess_init')]
    data[use_tbil1, -1] = data[use_tbil1, col_names.index('age_assess_repeat')]

    print("done")
    return data, col_names, ['bilirubin_age']


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


# TODO fix this for no df
'''
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
'''


def perform_association_subset(assoc_dir,
                               imputation_run_name,
                               region,
                               dep_var,
                               data,
                               col_names,
                               profile_mem_usage):
    # imputation_run_name = False means use imputed SNPs
    if profile_mem_usage:
        tracemalloc.start(50)
    chrom, poses = region.split(':')
    start, end = poses.split('-')
    region_string = f'{chrom}_{start}_{end}'
    with open(f'{assoc_dir}/run_logs/{dep_var}_{region_string}.log', 'w') as log, \
            open(f'{assoc_dir}/results/batches/{dep_var}_{region_string}.tab', 'w') as results:

        # if not profiling simply hide the function
        # behind a placeholder that does nothing
        if not profile_mem_usage:
            def tracemalloc_dump_snapshot(*args, **kwargs):
                pass
        else:
            def tracemalloc_dump_snapshot(fname, log, name):
                log.write(f'Memory usage at {name}: {tracemalloc.get_traced_memory()}\n')
                log.flush()
                #tracemalloc.take_snapshot().dump(fname)

        tracemalloc_dump_snapshot(f'{assoc_dir}/profiling/{region}_start', log, 'start')

        results.write("chrom\tpos\talleles\tlocus_info\tlocus_filtered\t#unrelated_samples_with_phenotype\t"
                      "#samples_GP_filtered\tallele_distribution\tmonoallelic_skipped"
                      f"\tp_{dep_var}\tcoeff_{dep_var}\tcoeff_intercept\n")
        results.flush()

        n_loci = 0
        batch_time = 0
        batch_size = 50
        total_time = 0

        use_strs = True
        if not imputation_run_name:
            use_strs = False

        if use_strs:
            genotype_iter = load_and_filter_genotypes.filtered_strs(
                imputation_run_name, region
            )
            # skip the samples, trust that they are in the same order
            next(genotype_iter)
        else:
            genotype_iter = load_and_filter_genotypes.filtered_imputed_snps(
                region
            )
            #samples should already be in the correct order
            #as it is prespecified to be that in the .sample file

        n_samples_with_phenotype = \
                np.sum(~np.isnan(data[:, col_names.index(f'{dep_var}_residual')]))

        start_time = time.time()
        tracemalloc_dump_snapshot(f'{assoc_dir}/profiling/{region}_pre_loop',
                         log, 'pre_loop')
        for len_gts, chrom, pos, allele_names, locus_info, locus_filtered in genotype_iter:
            tracemalloc_dump_snapshot(f'{assoc_dir}/profiling/{region}_loop_{pos}_start',
                            log, f'loop_{pos}_start')
            n_loci += 1
            results.write(f"{chrom}\t{pos}\t{allele_names}\t{locus_info}\t")
            if locus_filtered:
                results.write(f'{locus_filtered}\tnan\tnan\tnan\tnan\t1\tnan\tnan\n')
                results.flush()
                continue
            else:
                results.write('False\t')

            avg_len_gt = np.sum(len_gts, axis=1)/2
            gt = avg_len_gt

            unfiltered_samples = (~np.isnan(gt) &
                                  ~np.isnan(data[:, col_names.index(f'{dep_var}_residual')]))
            n_filtered_alleles = n_samples_with_phenotype - np.sum(unfiltered_samples)
            results.write(f"{n_samples_with_phenotype}\t{n_filtered_alleles}\t")

            alleles = np.unique(gt[unfiltered_samples], return_counts = True)
            any_alleles = False
            for count, allele in enumerate(alleles[0]):
                any_alleles = True
                if count>0:
                    results.write(",")
                results.write(f"{allele}:{alleles[1][count]}")
            if not any_alleles:
                results.write("nan")
            results.write("\t")

            if len(alleles[0]) <= 1:
                results.write("True\t1\tnan\tnan\n")
                results.flush()
                continue
            else:
                results.write("False\t")
                results.flush()
            tracemalloc_dump_snapshot(f'{assoc_dir}/profiling/{region}_loop_{pos}_preOLS',
                             log, f'loop_{pos}_preOLS')


            gt_const = np.concatenate((
                gt[unfiltered_samples].reshape(-1, 1),  np.ones((np.sum(unfiltered_samples), 1))
            ), axis = 1)

            #do da regression
            model = OLS(
                data[unfiltered_samples, col_names.index(f'{dep_var}_residual')],
                gt_const,
                missing='drop'
            )
            reg_result = model.fit()
            tracemalloc_dump_snapshot(f'{assoc_dir}/profiling/{region}_loop_{pos}_postOLS',
                             log, f'loop_{pos}_postOLS')
            pval = reg_result.pvalues[0]
            coef = reg_result.params[0]
            intercept_coef = reg_result.params[1]
            results.write(f"{pval:.2e}\t{coef}\t{intercept_coef}\n")
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
        if n_loci > 0:
            log.write(
                f"Done.\nTotal loci: {n_loci}\nTotal time: {total_time}s\ntime/locus: {total_time/n_loci}s\n"
            )
        else:
            log.write(f"No variants found in the region {region}\n")
        log.flush()
    return region_string


def invoke_plink(assoc_dir, chrom, dep_var, indep_vars):
    with open(f'{assoc_dir}/run_logs/{dep_var}_{chrom}.log', 'w') as log:
        plink_run_dir = f'{assoc_dir}/results/chr{chrom}'
        command = f"""
        {{
        cd {plink_run_dir}  || {{ echo "Failed to move to \
        {plink_run_dir}" ; exit 1 ; }} ; 
        {ukb}/utilities/plink2 \
            --pheno {assoc_dir}/covars_and_phenotypes.tab \
            --no-psam-pheno \
            --pheno-name {dep_var} \
            --covar-name {' '.join(indep_vars)} \
            --pfile {ukb}/array_imputed/pfile_converted/chr{chrom} \
            --chr {chrom} \
            --mac 20 \
            --glm omit-ref pheno-ids intercept \
            --ci 0.99999995 \
            --memory 128000 \
            --extract {ukb}association/temp/varnames.txt \
            --threads 28
        }} \
             > {assoc_dir}/run_logs/{dep_var}_{chrom}.plink.stdout \
            2> {assoc_dir}/run_logs/{dep_var}_{chrom}.plink.stderr
        """
        log.write('Command: ' + command)
        log.flush()
    try:
        sp.run(command, check=True, shell=True)
    except Exception as e:
        raise RuntimeError("plink failed!") from e


def run_plink_associations(assoc_dir, dep_vars, indep_vars, pheno_indep_vars, chromosome=None):
    dask_output_dir = f'{assoc_dir}/dask_output_logs'
    os.mkdir(dask_output_dir)
    os.mkdir(f'{assoc_dir}/run_logs')
    os.mkdir(f'{assoc_dir}/results')

    if chromosome:
        chroms = {int(chromosome)}
    else:
        chroms = range(1, 23)
    for chrom in chroms:
        os.mkdir(f'{assoc_dir}/results/chr{chrom}')

    cluster = dask_jobqueue.PBSCluster(
        queue="hotel",
        name="linear_associations",
        walltime="72:00:00",
        cores=28,
        memory="128GB",
        processes=1,
        log_directory=dask_output_dir
    )
    cluster.adapt(maximum_jobs=len(chroms)*len(dep_vars))
    client = dask.distributed.Client(cluster)

    print("Writing out results for each phenotype, chrom pair. "
          "see run_logs/<phenotype>_<chrom>.log for commands executed, "
          "run_logs/<phenotype>_<chrom>.plink.std<err/out> for plink output "
          "and see results/chr<chrom>/ for results.",
          flush=True)

    tasks = []
    for dep_var in dep_vars:
        curr_indep_vars = indep_vars.copy()
        curr_indep_vars.extend(pheno_indep_vars[dep_var])
        for chrom in chroms:
            tasks.append(dask.delayed(invoke_plink)(
                assoc_dir,
                chrom,
                dep_var,
                curr_indep_vars
            ))
    futures = client.compute(tasks)
    for future in futures:
        future.result() #block till all code is done executing
        # see here
        # https://distributed.dask.org/en/latest/manage-computation.html
    # TODO concat results across chroms here?


def concat_batches(results_dir, dep_var, region_strings):
    with open(f'{results_dir}/{dep_var}.tab', 'w') as outfile:
        first = True
        for string in region_strings:
            with open(f'{results_dir}/batches/{dep_var}_{string}.tab') as infile:
                first_line = True
                for line in infile:
                    if first_line:
                        first_line = False
                        if first:
                            first = False
                            outfile.write(line)
                        continue
                    outfile.write(line)

def run_associations_few_loci(assoc_dir, imputation_run_name, data, col_names, dep_vars,
                              loci, profile_mem_usage=False):
    os.mkdir(f'{assoc_dir}/run_logs')
    os.mkdir(f'{assoc_dir}/results')
    os.mkdir(f'{assoc_dir}/results/batches')
    if profile_mem_usage:
        os.mkdir(f'{assoc_dir}/profiling')
   
    print("Writing out results for each phenotype, region pair. "
          "see run_logs/<phenotype>_<region>.log for logs and "
          "see results/batches/<phenotype>_<region>.tab for results.",
          flush=True)
    for dep_var in dep_vars:
        locus_strings = []
        for locus in loci:
            perform_association_subset(
                assoc_dir,
                imputation_run_name,
                locus,
                dep_var,
                data,
                col_names,
                profile_mem_usage
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


def run_associations(assoc_dir, imputation_run_name, data, col_names, dep_vars,
                     chromosome=None, profile_mem_usage=False):
    dask_output_dir = f'{assoc_dir}/dask_output_logs'
    os.mkdir(dask_output_dir)
    os.mkdir(f'{assoc_dir}/run_logs')
    os.mkdir(f'{assoc_dir}/results')
    os.mkdir(f'{assoc_dir}/results/batches')
    if profile_mem_usage:
        os.mkdir(f'{assoc_dir}/profiling')
    cluster = dask_jobqueue.PBSCluster(
        queue="hotel",
        name="linear_associations",
        walltime="4:00:00",
        log_directory=dask_output_dir,
        processes=1,
        cores=2
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
        contents = f.read()
    client.register_worker_callbacks(
        setup=functools.partial(
            _worker_upload, data=contents, fname=fname,
        )
    )


    print("Writing out results for each phenotype, region pair. "
          "see run_logs/<phenotype>_<region>.log for logs and "
          "see results/batches/<phenotype>_<region>.tab for results.",
          flush=True)

    # prep the data frame for being distributed across the cluster
    # see
    # https://docs.dask.org/en/latest/delayed-best-practices.html#avoid-repeatedly-putting-large-inputs-into-delayed-calls
    data = dask.delayed(data)

    #get chrom lengths
    chr_lens = np.genfromtxt(
        f'{ukb}/misc_data/genome/chr_lens.txt',
        usecols=[1],
        skip_header=1,
        dtype=int
    )
    if chromosome:
        chroms = {int(chromosome)}
    else:
        chroms = range(1, 23)

    # many more imputed SNPs than imputed STRs
    # so to capture ~ the same number per task
    # use a smaller window for the former
    if imputation_run_name:
        # STR run
        range_per_task = int(1e7)
    else:
        # imputed SNP run
        range_per_task = int(1e6)

    overall_dep_tasks = []
    for dep_var in dep_vars:
        tasks = []
        for chrom in chroms:
            for start_pos in range(1, int(chr_lens[chrom-1]), range_per_task):
                tasks.append(dask.delayed(perform_association_subset)(
                    assoc_dir,
                    imputation_run_name,
                    f'{chrom}:{start_pos}-{start_pos + range_per_task - 1}',
                    dep_var,
                    data,
                    col_names,
                    profile_mem_usage
                ))
        overall_dep_tasks.append(dask.delayed(concat_batches)(
            f'{assoc_dir}/results', dep_var, tasks
        ))
    futures = client.compute(overall_dep_tasks)
    for future in futures:
        future.result() #block till all code is done executing
        # see here
        # https://distributed.dask.org/en/latest/manage-computation.html


def get_residuals_linear(readme,
                         assoc_dir,
                         data,
                         col_names,
                         dep_vars,
                         indep_vars,
                         pheno_indep_vars):
    """
    Linearlly regress out covariates

    col_names will be modified

    Returns
    -------
    data :
        To be used instead of the input data object
    """
    readme.write("Performing linear association against listed covariates"
                 " to get phenotypes residuals\n")
    readme.flush()

    print("Running linear associations of phenotypes against covariates to get"
          " residuals ... ", end='', flush=True)
    for dep_var in dep_vars:
        curr_indep_vars = indep_vars.copy()
        curr_indep_vars.extend(col_names.index(covar) for covar in pheno_indep_vars[dep_var])
        y = data[:, col_names.index(dep_var)]
        X = data[:, curr_indep_vars]
        with open(f'{assoc_dir}/{dep_var}_residual_assoc_summary.tab', 'w') as res_results:
            model = OLS(y, X, missing='drop')
            reg = model.fit()
            res_results.write("indep_var\tp\tcoeff\n")
            for idx, var in enumerate(curr_indep_vars):
                res_results.write(f'{var}\t{reg.pvalues[idx]:.2e}\t'
                                  f'{reg.params[idx]}\n')
        predictions = reg.predict(data[:, curr_indep_vars])
        readme.write("Training root mean square error for linear regression "
                     f"of phenotype {dep_var} against covariates: ")
        not_nan = ~np.isnan(y) & np.all(~np.isnan(X), axis=1)
        readme.write(str(
            np.std(data[not_nan, col_names.index(dep_var)] - predictions[not_nan])
        ))
        readme.write("\n")
        readme.flush()
        
        data = np.concatenate((
            data,
            np.full((data.shape[0], 1), np.nan)
        ), axis=1)
        col_names.append(f'{dep_var}_residual')
        data[:, -1] = y - predictions
    print("done", flush=True)
    return data


def get_residuals_csaps(readme,
                        assoc_dir,
                        data,
                        col_names,
                        dep_vars,
                        indep_vars,
                        pheno_indep_vars):
    """
    Regress out covariates with csaps

    col_names will be modified

    Returns
    -------
    data :
        To be used instead of the input data object
    """
    kfold = sklearn.model_selection.KFold(
        n_splits=10,
        shuffle=True,
        random_state=(hash('csaps') % 2**32)
    )
    rng = numpy.random.default_rng(np.abs(hash('csaps')))

    smooth_choices = np.array(
        [0, 1e-10, 3e-10, 1e-9, 3e-9,
            1e-8,  3e-8,  1e-7, 3e-7,
            1e-6,  3e-6,  1e-5, 3e-5,
            1e-4,  3e-4,  1e-3, 3e-3,
            1e-2,  3e-2,  1e-1, 3e-1])
    """
    smooth_choices = np.array([0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 3e-1])
    """
    smooth_choices = list(smooth_choices) + [0.5] + list(1 - smooth_choices)

    readme.write("Performing cubic smoothing spline associations against listed"
                 " covariates to get phenotypes residuals\n")
    readme.write("Starting with residuals as phenotype values centered, per sex, with mean zero. "
                 "Sequentially taking each covariate in turn, tuning its smoothing parameter "
                 "using 5-fold cross validation, and upon choosing the parameter value recalculating "
                 "the residuals by regressing out that covariate. (Parameter ranges from 0 (least square lines) "
                 "to 1 (perfect fit), where "
                 "increasing values allow for cubics with greater curavture). This model assumes that "
                 "contributions from each covariate are additive. This model requires that covariate values be "
                 "unique, so a small jitter (magnitude <= 5e-4) is added to each covaraite value during model fitting. "
                 "Smoothing parameters are being chosen from: ")
    readme.write(str(smooth_choices))
    readme.write('\n')
    readme.flush()

    print("Running CSAP associations of phenotypes against covariates to get"
          " residuals ... ", end='', flush=True)
    for dep_var in dep_vars:
        curr_indep_vars = indep_vars.copy()
        curr_indep_vars.extend(col_names.index(covar) for covar in pheno_indep_vars[dep_var])
        curr_indep_vars.remove(col_names.index('intercept'))
        curr_indep_vars.remove(col_names.index("sex"))
        smooth_param = [None]*len(curr_indep_vars)

        readme.write(f"Covariates in order for dep var {dep_var}: ")
        readme.write(str([col_names[idx] for idx in curr_indep_vars]))
        readme.write('\n')
        readme.flush()

        y = data[:, col_names.index(dep_var)]
        X = data[:, curr_indep_vars]
        not_nan = ~np.isnan(y) & np.all(~np.isnan(X), axis=1)
        y = y[not_nan]
        X = X[not_nan, :]
        # jitter any lower than this causes some values to end up being
        # duplicates due to loss of precision in 64bit floats which causes
        # CSAP to crash
        jitter = 5e-4*(rng.random(X.shape) - 0.5)
        jittered_X = X + jitter
        curr_residuals = y.copy()

        # starting residuals
        sex_idx = col_names.index("sex")
        sex_is_1 = data[not_nan, sex_idx] == 1
        curr_residuals[sex_is_1] -= np.mean(curr_residuals[sex_is_1])
        curr_residuals[~sex_is_1] -= np.mean(curr_residuals[~sex_is_1])
        curr_rmse = np.std(curr_residuals)
        print("Starting residual rmse: ", curr_rmse)

        with open(f'{assoc_dir}/{dep_var}_residual_assoc_summary.txt', 'w') as res_results:
            res_results.write("Covariate\tSmoothing Parameter\n")
            for idx in range(len(curr_indep_vars)):
                best_rmse = curr_rmse
                best_smooth_choice = None
                for smooth_choice in smooth_choices:
                    rmse = 0
                    for train, test in kfold.split(jittered_X):
                        sort = np.argsort(jittered_X[train, idx])
                        spline = csaps.csaps(jittered_X[train, idx][sort],
                                             curr_residuals[train][sort],
                                             smooth=smooth_choice)
                        print(np.std(curr_residuals[test] - spline(X[test, idx])))
                        rmse += np.std(curr_residuals[test] - spline(X[test, idx]))
                    rmse /= kfold.get_n_splits()
                    print(f"Avg rmse for smoothing param {smooth_choice}: {rmse}")
                    if rmse <= best_rmse:
                        best_rmse = rmse
                        best_smooth_choice = smooth_choice
                smooth_param[idx] = best_smooth_choice
                print(f"Best avg. test rmse for idx {idx}: {best_rmse}")
                print(f"Smoothing param: {best_smooth_choice}")
                # update current residuals by fitting to the entire dataset
                if best_smooth_choice:
                    sort = np.argsort(jittered_X[:, idx])
                    spline = csaps.csaps(jittered_X[sort, idx],
                                         curr_residuals[sort],
                                         smooth=best_smooth_choice)
                    predictions = spline(X[:, idx])
                    curr_residuals = curr_residuals - predictions
                    curr_rmse = np.std(curr_residuals)
                    print(f"Current RMSE after regression on idx {idx}: ", curr_rmse)
                else:
                    print(f"Not regression on idx {idx}")
        # now curr_residuals are the final residuals and curr_rmse
        # is the final rmse
        readme.write(f"Last validation root mean square error for phenotype {dep_var} "
                     "smoothing spline regression against covariates: ")
        readme.write(str(best_rmse))
        readme.write("\n")
        readme.write(f"Smoothing parameters for {dep_var}: ")
        readme.write(str(smooth_param))
        readme.write("\n")
        readme.write("(Smoothing parameters with value None indicate that the covariate was not used during prediction)\n")
        readme.write(f"Training root mean square error for phenotype {dep_var} "
                     "smoothing spline regression against covariates: ")
        readme.write(str(curr_rmse))
        readme.write('\n')
        readme.flush()

        data = np.concatenate((
            data,
            np.full((data.shape[0], 1), np.nan)
        ), axis=1)
        col_names.append(f'{dep_var}_residual')
        data[not_nan, -1] = curr_residuals
    print("done", flush=True)
    return data

def main():  # noqa: D103
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_filtering_run_name')
    parser.add_argument(
        'association_run_name',
        help="The name to give this association run"
    )
    parser.add_argument('--imputation_run_name')
    parser.add_argument(
        '--loci',
        help='comma separated list with elements of the form chr:pos'
    )
    parser.add_argument(
        '--debug',
        help='omits some steps to be faster',
        action='store_true'
    )
    parser.add_argument(
        '--no-run',
        help="do everything up to running the association, but don't run it",
        action='store_true'
    )
    parser.add_argument(
        '--imputed-snps',
        action='store_true',
        help="Use imputed SNP genotypes instead of imputed STR genotypes"
    )
    parser.add_argument(
        '--plink',
        action="store_true",
        help=("Use plink for computing SNP associations. Only meaningful "
              "if --imputed-snps is set")
    )
    parser.add_argument(
        '--chromosome',
        help='run linear associations only for one chrom instead of genome wide'
    )
    parser.add_argument(
        '--profile-memory-usage',
        action='store_true'
    )
    parser.add_argument(
        '--csaps',
        help='use csaps for covariate regression instead of linear regression',
        action='store_true'
    )
    args = parser.parse_args()

    assert not (args.chromosome and args.loci)

    if args.plink:
        assert args.imputed_snps

    if not (
        (args.imputation_run_name and not args.imputed_snps) or
        (not args.imputation_run_name and args.imputed_snps)
    ):
        print('Error: Expecting exactly one of the arguments --imputation-run-name '
              'and --imputed-snps')
        sys.exit()

    sample_filtering_run_name = args.sample_filtering_run_name
    samp_fname = f'{ukb}/sample_qc/runs/{sample_filtering_run_name}/combined_unrelated.sample'
    if not os.path.exists(samp_fname):
        print(f"Sample file {samp_fname} does not exist!")
        sys.exit(1)

    use_strs = not args.imputed_snps
    if use_strs:
        imputation_run_name = args.imputation_run_name
        assert os.path.exists(f'{ukb}/str_imputed/runs/{imputation_run_name}/')

    association_run_name = args.association_run_name
    assoc_dir = f'{ukb}/association/runs/{association_run_name}'
    if os.path.exists(assoc_dir):
        print(f"Association with run name {association_run_name} already "
              f"exists!", file=sys.stderr)
        sys.exit(1)

    os.mkdir(assoc_dir)

    if args.loci:
        pre_format_loci = args.loci.split(',')
        loci = []
        for locus in pre_format_loci:
            chrom, pos = locus.split(':')
            loci.append(f'{chrom}:{pos}-{pos}')
    else:
        loci = None


    with open(f"{assoc_dir}/README", 'w') as readme, \
            open(f"{assoc_dir}/phenotype_units.txt", 'w') as phenotypes:

        readme.write(f"Run args: {args}\n")
        today = datetime.datetime.now().strftime("%Y_%m_%d")
        readme.write(f"Run date: {today}\n")
        readme.flush()
        phenotypes.write("phenotype:unit\n")
        phenotypes.flush()

        data, col_names, indep_vars = load_covars(readme)
        data, col_names, height_indep_vars = load_height(data, col_names, readme, phenotypes)
        data, col_names, bilirubin_indep_vars = load_bilirubin(data, col_names, readme, phenotypes)

        # covariate columns:
        dep_vars = [
            'height',
            'total_bilirubin'
        ]
        pheno_indep_vars = {
            'height': height_indep_vars,
            'total_bilirubin': bilirubin_indep_vars
        }

        print("Setting up samples ... ", end="", flush=True)
        # Already calculated the largest unrelated subset of a filtering list
        # based on height. Will need to do this more intelligently when we
        # don't use the same samples for bilirubin and height
        floc = f'{ukb}/sample_qc/runs/{sample_filtering_run_name}/combined_unrelated.sample'
        readme.write(f"Subsetting to samples from file {floc}.\n")
        readme.flush()
        sample_subset = np.genfromtxt(
            floc,
            usecols=[0],
            skip_header=1,
            delimiter=" "
        ).reshape(-1, 1)
        data = merge_arrays(sample_subset, data)

        # order the observations in data
        # based on the order of samples in chr1 vcf
        if use_strs:
            vcf = cyvcf2.VCF(
                f'{ukb}/str_imputed/runs/{imputation_run_name}/'
                f'vcfs/annotated_strs/chr1.vcf.gz'
            )
            samples_array = np.array(
                list(int(sample.split("_")[0]) for sample in vcf.samples)
            )
        else:
            bgen_samples = []
            with open(f'{ukb}/microarray/ukb46122_hap_chr1_v2_s487314.sample') as samplefile:
                for num, line in enumerate(samplefile):
                    if num <= 1:
                        # skip first two lines
                        continue
                    bgen_samples.append(line.split()[0])
            assert len(bgen_samples) == 487409
            samples_array = np.array(bgen_samples, dtype=float)
        samples_array = samples_array.reshape(-1, 1)
        data = merge_arrays(samples_array, data)
        print("done")

        # Calculate the number of samples used for each phenotype
        for phen in dep_vars:
            readme.write(f'Number of samples for phenotype {phen}: ')
            nonnans = np.all(~np.isnan(data[:, indep_vars]), axis=1)
            nonnans &= ~np.isnan(data[:, col_names.index(phen)])
            for covar in pheno_indep_vars[phen]:
                nonnans &= ~np.isnan(data[:, col_names.index(covar)])
            readme.write(str(np.sum(nonnans)))
            readme.write('\n')
            readme.flush()

        if not args.debug:
            plot_dirs = {}
            for dep_var in dep_vars:
                plot_dirs[dep_var] = f'{assoc_dir}/{dep_var}_summary_plots'
                os.mkdir(plot_dirs[dep_var])

        if not args.debug:
            plot_phenotype_by_sex(
                df,
                'total_bilirubin',
                f'{plot_dirs["total_bilirubin"]}/measured.png',
                'umol/L',
                'Measured total_bilirubin'
            )

        readme.write("Log transforming total_bilirubin\n")
        readme.flush()
        # figure out how to tie this to the unit it the units file above
        data[:, col_names.index('total_bilirubin')] = \
                np.log(data[:, col_names.index('total_bilirubin')])

        if not args.debug:
            plot_phenotype_by_sex(
                data,
                'total_bilirubin',
                f'{plot_dirs["total_bilirubin"]}/log_transform.png',
                'log(umol/L)',
                'Log transform total_bilirubin'
            )

            plot_phenotype_by_sex(
                data,
                'height',
                f'{plot_dirs["height"]}/measured.png',
                'cm',
                'Measured height'
            )

        units = {'height': 'cm', 'total_bilirubin': 'log(umol/L)'}

        for dep_var in dep_vars:
            readme.write(f"Standard deviation of {dep_var}: " )
            idx = col_names.index(dep_var)
            not_nan = ~np.isnan(data[:, idx])
            for cov_idx in indep_vars:
                not_nan &= ~np.isnan(data[:, cov_idx])
            for cov in pheno_indep_vars[dep_var]:
                not_nan &= ~np.isnan(data[:, col_names.index(cov)])
            readme.write(str(np.std(data[not_nan, idx])))
            readme.write("\n")

        if not args.plink:
            if not args.csaps:
                data = get_residuals_linear(
                    readme,
                    assoc_dir,
                    data,
                    col_names,
                    dep_vars,
                    indep_vars,
                    pheno_indep_vars
                )
            else:
                data = get_residuals_csaps(
                    readme,
                    assoc_dir,
                    data,
                    col_names,
                    dep_vars,
                    indep_vars,
                    pheno_indep_vars
                )

        assert len(col_names) == data.shape[1]

        print("Writing out sample x {phenotypes, covaraites} tab file ... ",
              end="", flush=True)
        #Write out the csv we're going to be running associations on
        # use plink syntax for ID and FID in case we want to run plink associations
        dcopy = np.concatenate((data[:, 0:1], data), axis=1)
        cols_copy = col_names.copy()
        cols_copy.insert(0, 'FID')
        cols_copy[1] = 'IID'
        np.savetxt(
            f'{assoc_dir}/covars_and_phenotypes.tab',
            dcopy,
            delimiter='\t',
            header='\t'.join(cols_copy),
            fmt='%.16g',
            comments='#'
        )
        del dcopy
        print("done", flush=True)

        if not args.debug:
            for dep_var in dep_vars:
                plot_phenotype_by_sex(
                    data,
                    f'{dep_var}_residual',
                    f'{plot_dirs[dep_var]}/residual.png',
                    f'residual {units[dep_var]}',
                    f'Residual {dep_var}'
                )
        if args.no_run:
            readme.write('--no-run specified. Exiting\n')
            readme.flush()
            sys.exit()
        if use_strs:
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
                         "increase of 1 repeat in both haplotypes.\n")
            readme.flush()

            if loci is not None:
                run_associations_few_loci(assoc_dir, imputation_run_name,
                                          data, col_names, dep_vars, loci,
                                          profile_mem_usage=args.profile_memory_usage)
            else:
                run_associations(assoc_dir, imputation_run_name, data,
                                 col_names, dep_vars,
                                chromosome=args.chromosome,
                                profile_mem_usage=args.profile_memory_usage)
        elif not args.plink:
            readme.write("Using imputed snp calls\n")
            readme.write(
                "Should filter loci whose INFO scores are <0.3 "
                "https://www.ukbiobank.ac.uk/wp-content/uploads/2014/04/imputation_documentation_May2015.pdf"
                " (or <0.8 "
                "http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas"
                " ?). Not doing this yet to let it be tunable downstream.\n"
            )
            readme.write(
                "Filtering calls whose probability as estimated by the UKB "
                "imputation algorithm is less than 0.9. Number filtered this "
                "way is the #samples_GP_filtered column in "
                "filtering/<region>.txt\n"
            )
            readme.write(
                "Performing linear association of residual phenotypes "
                "against averaged genotype (0 homozygous ref, 0.5 "
                "heterozygous, 1 homozygous alt). "
                "Coefficients are reported as change in "
                "phenotype when changing both haplotypes from the 0 "
                "allele to the 1 allele\n"
            )
            readme.flush()
            if loci is not None:
                run_associations_few_loci(assoc_dir, False, data, col_names, dep_vars, loci,
                                          profile_mem_usage=args.profile_memory_usage)
            else:
                run_associations(assoc_dir, False, data, col_names, dep_vars,
                                 chromosome=args.chromosome,
                                 profile_mem_usage=args.profile_memory_usage)
        else:
            readme.write(
                "Using imputed snp calls\n"
                "Running associations with plink\n"
                "Need to filter loci with e.g. <20 nonmajor alleles. Recommended by plink docs "
                "https://www.cog-genomics.org/plink/2.0/assoc and by ukb "
                "GWAS summary "
                "http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas"
                " but have not done so yet so this choice can be tuned downstream.\n"
                "Genotypes are encoded as 0 - ref, 1 - alt\n"
                "Running linear assocation with all relevant covariates "
                "included (no regressing out the covariates first and then "
                "comparing the genotype to the phenotype)\n"
            )
            readme.flush()
            indep_vars.remove('intercept')
            run_plink_associations(assoc_dir, dep_vars, indep_vars, pheno_indep_vars,
                                   chromosome=args.chromosome)

if __name__ == "__main__":
    main()

