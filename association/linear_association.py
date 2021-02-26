
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
import traceback
import tracemalloc
from typing import List, Optional

import cyvcf2
import dask
import dask.distributed
import dask_jobqueue
import matplotlib.pyplot as plt
import numpy as np
import numpy.random
import scipy.stats
import sklearn.ensemble
import sklearn.kernel_ridge
import sklearn.model_selection
import sklearn.neighbors
from statsmodels.regression.linear_model import OLS
import statsmodels.stats.weightstats

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
    readme.write(f"{ids_and_sex.shape[0]} total participants")
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
        try:
            perform_association_subset_helper(
                assoc_dir,
                imputation_run_name,
                region,
                dep_var,
                data,
                col_names,
                profile_mem_usage,
                log,
                results
            )
        except Exception as e:
            log.write(traceback.format_exc())
            log.flush()
            raise e
    return region_string


def perform_association_subset_helper(
        assoc_dir,
        imputation_run_name,
        region,
        dep_var,
        data,
        col_names,
        profile_mem_usage,
        log,
        results):
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

    results.write("chrom\tpos\talleles\tlocus_filtered\t"
                  f"p_{dep_var}\tcoeff_{dep_var}\tcoeff_intercept\t")
    results.flush()

    n_loci = 0
    batch_time = 0
    batch_size = 50
    total_time = 0

    use_strs = True
    if not imputation_run_name:
        use_strs = False

    unfiltered_samples = ~np.isnan(data[:, col_names.index(f'{dep_var}_inv_norm_rank_residual')])

    if use_strs:
        genotype_iter = load_and_filter_genotypes.load_strs(
            imputation_run_name, region, unfiltered_samples
        )
    else:
        genotype_iter = load_and_filter_genotypes.load_imputed_snps(
            region, unfiltered_samples
        )

    # first yield is special
    extra_detail_fields = next(genotype_iter)
    results.write('\t'.join(extra_detail_fields) + '\t')
    results.write(f'mean_residual_{dep_var}_per_single_dosage\t'
                  '0.05_significance_CI\t'
                  '5e-8_significance_CI')

    if use_strs:
        results.write(f'\tmean_residual_{dep_var}_per_paired_dosage\t'
                      '0.05_significance_CI\t'
                      '5e-8_significance_CI')
    results.write('\n')
    results.flush()

    start_time = time.time()
    tracemalloc_dump_snapshot(f'{assoc_dir}/profiling/{region}_pre_loop',
                     log, 'pre_loop')
    for dosage_gts, unique_alleles, chrom, pos, locus_filtered, locus_details in genotype_iter:
        assert len(locus_details) == len(extra_detail_fields)
        tracemalloc_dump_snapshot(f'{assoc_dir}/profiling/{region}_loop_{pos}_start',
                        log, f'loop_{pos}_start')

        n_loci += 1
        allele_names = ','.join(list(unique_alleles.astype(str)))
        results.write(f"{chrom}\t{pos}\t{allele_names}\t")
        if locus_filtered:
            results.write(f'{locus_filtered}\t1\tnan\tnan\t')
            results.write('\t'.join(locus_details))
            if use_strs:
                results.write('\tnan'*6 + '\n')
            else:
                results.write('\tnan'*3 + '\n')
            results.flush()
            continue
        else:
            results.write('False\t')

        tracemalloc_dump_snapshot(f'{assoc_dir}/profiling/{region}_loop_{pos}_preOLS',
                         log, f'loop_{pos}_preOLS')

        if use_strs:
            gts = np.sum([_len*np.sum(dosages, axis=1) for
                          _len, dosages in dosage_gts.items()], axis=0)
        else:
            gts = dosage_gts[:, 1] + 2*dosage_gts[:, 2]
        gt_const = np.concatenate((
            gts.reshape(-1, 1),  np.ones((gts.shape[0], 1))
        ), axis = 1)

        #do da regression
        model = OLS(
            data[unfiltered_samples, col_names.index(f'{dep_var}_inv_norm_rank_residual')],
            gt_const,
            missing='drop'
        )
        reg_result = model.fit()
        tracemalloc_dump_snapshot(f'{assoc_dir}/profiling/{region}_loop_{pos}_postOLS',
                         log, f'loop_{pos}_postOLS')
        pval = reg_result.pvalues[0]
        coef = reg_result.params[0]
        intercept_coef = reg_result.params[1]
        results.write(f"{pval:.2e}\t{coef}\t{intercept_coef}\t")
        results.write('\t'.join(locus_details) + '\t')

        if use_strs:
            single_dosages = {}
            paired_dosages = {}
            for len1 in unique_alleles:
                for len2 in unique_alleles:
                    if len1 != len2:
                        dosages = (dosage_gts[len1][:, 0]*dosage_gts[len2][:, 1] +
                                   dosage_gts[len1][:, 1]*dosage_gts[len2][:, 0])
                    else:
                        dosages = dosage_gts[len1][:, 0]*dosage_gts[len1][:, 1]
                    if len1 + len2 not in single_dosages:
                        single_dosages[len1 + len2] = dosages
                    else:
                        single_dosages[len1 + len2] += dosages
                    minlen = min(len1, len2)
                    maxlen = max(len1, len2)
                    paired_dosages[(minlen, maxlen)] = dosages
            single_dosage_means = {}
            single_dosage_95_CI = {}
            single_dosage_GWAS_CI= {}
            for _len, dosages in single_dosages.items():
                mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                    data[unfiltered_samples, col_names.index(f'{dep_var}_residual')],
                    weights = dosages
                )
                single_dosage_means[_len] = mean_stats.mean
                single_dosage_95_CI[_len] = mean_stats.tconfint_mean()
                single_dosage_GWAS_CI[_len] = mean_stats.tconfint_mean(5e-8)
            paired_dosage_means = {}
            paired_dosage_95_CI = {}
            paired_dosage_GWAS_CI= {}
            for _len, dosages in paired_dosages.items():
                mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                    data[unfiltered_samples, col_names.index(f'{dep_var}_residual')],
                    weights = dosages
                )
                paired_dosage_means[_len] = mean_stats.mean
                paired_dosage_95_CI[_len] = mean_stats.tconfint_mean()
                paired_dosage_GWAS_CI[_len] = mean_stats.tconfint_mean(5e-8)
            results.write(load_and_filter_genotypes.dict_str(single_dosage_means) + '\t')
            results.write(load_and_filter_genotypes.dict_str(single_dosage_95_CI) + '\t')
            results.write(load_and_filter_genotypes.dict_str(single_dosage_GWAS_CI) + '\t')
            results.write(load_and_filter_genotypes.dict_str(paired_dosage_means) + '\t')
            results.write(load_and_filter_genotypes.dict_str(paired_dosage_95_CI) + '\t')
            results.write(load_and_filter_genotypes.dict_str(paired_dosage_GWAS_CI) + '\n')
        else:
            single_dosage_means = {}
            single_dosage_95_CI = {}
            single_dosage_GWAS_CI= {}
            for alt_count in range(3):
                mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                    data[unfiltered_samples, f'{dep_var}_residual'],
                    weights = dosage_gts[:, alt_count]
                )
                single_dosage_means[alt_count] = mean_stats.mean
                single_dosage_95_CI[alt_count] = mean_stats.tconfint_mean()
                single_dosage_GWAS_CI[alt_count] = mean_stats.tconfint_mean(5e-8)
            results.write(load_and_filter_genotypes.dict_str(single_dosage_means) + '\t')
            results.write(load_and_filter_genotypes.dict_str(single_dosage_95_CI) + '\t')
            results.write(load_and_filter_genotypes.dict_str(single_dosage_GWAS_CI) + '\n')

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
            --pheno-name {dep_var}_inv_norm_rank \
            --covar-name {' '.join(indep_vars)} \
            --pfile {ukb}/array_imputed/pfile_converted/chr{chrom} \
            --chr {chrom} \
            --mac 20 \
            --glm omit-ref pheno-ids intercept \
            --ci 0.99999995 \
            --memory 128000 \
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


def run_plink_associations(assoc_dir, dep_vars, indep_vars, pheno_indep_vars, col_names, chromosome=None):
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
        curr_indep_vars = pheno_indep_vars[dep_var].copy()
        curr_indep_vars.extend([col_names[idx] for idx in indep_vars])
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

    train_size = 0.9
    splitter = sklearn.model_selection.ShuffleSplit(
        n_splits = 5,
        train_size = train_size,
        random_state = hash('linear') % 2**32,
    )

    print("Running linear associations of phenotypes against covariates to get"
          " residuals ... ", end='', flush=True)
    for dep_var in dep_vars:
        curr_indep_vars = indep_vars.copy()
        curr_indep_vars.extend(col_names.index(covar) for covar in pheno_indep_vars[dep_var])
        y = data[:, col_names.index(dep_var)]
        X = data[:, curr_indep_vars]
        not_nan = ~np.isnan(y) & np.all(~np.isnan(X), axis=1)
        y = y[not_nan]
        X = X[not_nan, :].copy()

        # estimate RMSE from linear model
        rmse = 0
        for train, test in splitter.split(X):
            model = OLS(y[train], X[train, :])
            reg = model.fit()
            rmse += np.sqrt(np.mean((y[test] - reg.predict(X[test]))**2))
        rmse /= splitter.get_n_splits()
        n_samples = y.shape[0]
        readme.write(f"Validation root mean square error for linear regression "
                     f"of phenotype {dep_var} "
                     f"using {train_size*n_samples:.0f} training samples and "
                     f"{(1-train_size)*n_samples:.0f} validation samples: {rmse:.4f}\n")

        with open(f'{assoc_dir}/{dep_var}_residual_assoc_summary.tab', 'w') as res_results:
            model = OLS(y, X)
            reg = model.fit()
            res_results.write("indep_var\tp\tcoeff\n")
            for idx, var in enumerate(curr_indep_vars):
                res_results.write(f'{var}\t{reg.pvalues[idx]:.2e}\t'
                                  f'{reg.params[idx]}\n')
        predictions = reg.predict(data[np.ix_(not_nan, curr_indep_vars)])

        data = np.concatenate((
            data,
            np.full((data.shape[0], 1), np.nan)
        ), axis=1)
        col_names.append(f'{dep_var}_residual')
        data[not_nan, -1] = y - predictions
        readme.flush()

    print("done", flush=True)
    return data


def get_residuals_adaboost(
        readme,
        assoc_dir,
        data,
        col_names,
        dep_vars,
        indep_vars,
        pheno_indep_vars):
    """
    Regress out covariates with random forests
    https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html#sklearn.ensemble.RandomForestRegressor

    col_names will be modified

    Returns
    -------
    data :
        To be used instead of the input data object
    """
    train_size = 0.9
    splitter = sklearn.model_selection.ShuffleSplit(
        n_splits=5,
        train_size=train_size,
        random_state=(hash('adaboost') % 2**32)
    )

    print("Running AdaBoost regression on phenotypes against covariates to get"
          " residuals ... ", flush=True)
    for dep_var in dep_vars:
        curr_indep_vars = indep_vars.copy()
        curr_indep_vars.extend(col_names.index(covar) for covar in pheno_indep_vars[dep_var])
        curr_indep_vars.remove(col_names.index('intercept'))

        y = data[:, col_names.index(dep_var)]
        X = data[:, curr_indep_vars]
        not_nan = ~np.isnan(y) & np.all(~np.isnan(X), axis=1)
        y = y[not_nan]
        X = X[not_nan, :].copy()
        n_samples = y.shape[0]

        print('Starting')
        n_stumps = 25
        count = 0
        while True:
            n_stumps *= 2

            rmse = 0
            start_time = time.time()
            for train, test in splitter.split(X):
                count += 1
                model = sklearn.ensemble.AdaBoostRegressor(
                    random_state = (hash('adaboost model')*count % 2**32),
                    n_estimators = n_stumps
                )
                model.fit(X[train, :], y[train])
                predictions = model.predict(X[test, :])
                rmse += np.sqrt(np.mean((y[test] - predictions)**2))
            rmse /= splitter.get_n_splits()
            readme.write(f"Validation root mean square error for phenotype {dep_var} "
                         f"AdaBoost regression using {n_stumps} stumps (n_estimators), "
                         f"{train_size*n_samples:.0f} training samples and "
                         f"{(1-train_size)*n_samples:.0f} validation samples: {rmse:.4f}\n")
            readme.flush()
            print(f"Time for building random forest model for {n_stumps} stumps: "
                  f"{time.time()-start_time:.0f}", flush=True)

        predictions = model.predict(X)
        residuals = y - predictions

        data = np.concatenate((
            data,
            np.full((data.shape[0], 1), np.nan)
        ), axis=1)
        col_names.append(f'{dep_var}_residual')
        data[not_nan, -1] = residuals
    print("done", flush=True)
    return data



def get_residuals_random_forest(
        readme,
        assoc_dir,
        data,
        col_names,
        dep_vars,
        indep_vars,
        pheno_indep_vars):
    """
    Regress out covariates with random forests
    https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html#sklearn.ensemble.RandomForestRegressor

    col_names will be modified

    Returns
    -------
    data :
        To be used instead of the input data object
    """
    train_size = 0.9
    splitter = sklearn.model_selection.ShuffleSplit(
        n_splits=5,
        train_size=train_size,
        random_state=(hash('random_forest') % 2**32)
    )

    print("Running random forest regression on phenotypes against covariates to get"
          " residuals ... ", flush=True)
    for dep_var in dep_vars:
        readme.write("Using random forest regression model with 10 min_samples_leaf\n")
        curr_indep_vars = indep_vars.copy()
        curr_indep_vars.extend(col_names.index(covar) for covar in pheno_indep_vars[dep_var])
        curr_indep_vars.remove(col_names.index('intercept'))

        y = data[:, col_names.index(dep_var)]
        X = data[:, curr_indep_vars]
        not_nan = ~np.isnan(y) & np.all(~np.isnan(X), axis=1)
        y = y[not_nan]
        X = X[not_nan, :].copy()
        n_samples = y.shape[0]

        start_time = time.time()
        n_trees = 50
        count = 0
        while True:
            n_trees *= 2

            rmse = 0
            for train, test in splitter.split(X):
                count += 1
                model = sklearn.ensemble.RandomForestRegressor(
                    min_samples_leaf = 10,
                    random_state = (hash('random_forest model')*count % 2**32),
                    n_jobs = 20,
                    n_estimators = n_trees
                )
                model.fit(X[train, :], y[train])
                predictions = model.predict(X[test, :])
                rmse += np.sqrt(np.mean((y[test] - predictions)**2))
            rmse /= splitter.get_n_splits()
            readme.write(f"Validation root mean square error for phenotype {dep_var} "
                         f"random forest regression using {n_trees} trees (n_estimators), "
                         f"{train_size*n_samples:.0f} training samples and "
                         f"{(1-train_size)*n_samples:.0f} validation samples: {rmse:.4f}\n")
            readme.flush()
            print(f"Time for building random forest model for {n_trees} trees: "
                  f"{time.time()-start_time:.0f}", flush=True)

        predictions = model.predict(X)
        residuals = y - predictions

        data = np.concatenate((
            data,
            np.full((data.shape[0], 1), np.nan)
        ), axis=1)
        col_names.append(f'{dep_var}_residual')
        data[not_nan, -1] = residuals
    print("done", flush=True)
    return data


def get_residuals_kernel_ridge(
        readme,
        assoc_dir,
        data,
        col_names,
        dep_vars,
        indep_vars,
        pheno_indep_vars):
    """
    Regress out covariates with Kernel Ridge Regression
    https://scikit-learn.org/stable/modules/generated/sklearn.kernel_ridge.KernelRidge.html#sklearn.kernel_ridge.KernelRidge
    (a slightly modified version of Support Vector Machine regression)

    col_names will be modified

    Returns
    -------
    data :
        To be used instead of the input data object
    """
    kfold = sklearn.model_selection.KFold(
        shuffle=True,
        random_state=(hash('krr') % 2**32)
    )
    rng = numpy.random.default_rng(np.abs(hash('krr')))
    subset_size = 50000 #data is too large to use all at once in training, so subset

    param_choices = [10**(i/8) for i in range(-80, -40)]

    readme.write(f"Training sample size for kernel ridge regression: {subset_size}\n")

    """
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
    """

    print("Running KRR on phenotypes against covariates to get"
          " residuals ... ", end='', flush=True)
    for dep_var in dep_vars:
        curr_indep_vars = indep_vars.copy()
        curr_indep_vars.extend(col_names.index(covar) for covar in pheno_indep_vars[dep_var])
        curr_indep_vars.remove(col_names.index('intercept'))
        print([col_names[idx] for idx in curr_indep_vars])

        y = data[:, col_names.index(dep_var)]
        X = data[:, curr_indep_vars]
        not_nan = ~np.isnan(y) & np.all(~np.isnan(X), axis=1)
        y = y[not_nan]
        X = X[not_nan, :].copy()
        # standard normalize the covariates
        X = (X - np.mean(X, axis=0))/np.std(X, axis=0)

        subset_idx = rng.choice(X.shape[0], size=subset_size, replace=False)
        y_subset = y[subset_idx]
        X_subset = X[subset_idx, :]

        best_rmse = np.inf
        best_params = (None, None)
        n_rounds = 0
        start_time = time.time()
        print('Starting', end='\r')
        for alpha in param_choices:
            for gamma in param_choices:
                rmse = 0
                for train, test in kfold.split(X_subset):
                    n_rounds += 1
                    model = sklearn.kernel_ridge.KernelRidge(alpha=alpha, kernel='rbf', gamma=gamma)
                    model.fit(X_subset[train, :], y_subset[train])
                    predictions = model.predict(X_subset[test, :])
                    rmse += np.sqrt(np.mean((y_subset[test] - predictions)**2))
                    curr_time = time.time()
                    print(f"\033[2K\033[1GModel fits so far: {n_rounds}. Avg time per fit: {(curr_time - start_time)/n_rounds}sec", end='\r', flush=True)
                rmse /= kfold.get_n_splits()
                print(f"\033[2K\033[1GAvg validation rmse for params (alpha, gamma) ({alpha:>7.4g}, {gamma:>7.4g}): {rmse:>7.4g}", flush=True)
                if rmse <= best_rmse:
                    best_rmse = rmse
                    best_params = (alpha, gamma)

        readme.write(f"Last validation root mean square error for phenotype {dep_var} "
                     "kernel ridge regression against covariates: ")
        readme.write(str(best_rmse))
        readme.write("\n")
        readme.write(f"Best KRR metaparameters for {dep_var}: ")
        readme.write(str(best_params))
        readme.write("\n")

        model = sklearn.kernel_ridge.KernelRidge(alpha=best_params[0],
                                                 kernel='rbf',
                                                 gamma=best_params[1])
        model.fit(X_subset, y_subset)
        predictions = model.predict(X)
        residuals = y - predictions
        rmse = np.sqrt(np.mean((residuals)**2))
        readme.write(f"Training root mean square error for phenotype {dep_var} "
                     "kernel ridge regression against covariates: ")
        readme.write(str(rmse))
        readme.write('\n')
        readme.flush()

        data = np.concatenate((
            data,
            np.full((data.shape[0], 1), np.nan)
        ), axis=1)
        col_names.append(f'{dep_var}_residual')
        data[not_nan, -1] = residuals
    print("done", flush=True)
    return data


def rank_phenotypes(
        readme,
        data,
        col_names,
        dep_vars,
        residuals=True):
    """
    col_names will be modified

    Ranks begin with 0
    """
    readme.write("Ranking phenotype residuals. Tie breaks for equal ranks are arbitrary and stable.\n")
    for dep_var in dep_vars:
        data = np.concatenate((
            data,
            np.full((data.shape[0], 1), np.nan)
        ), axis=1)

        if residuals:
            col_names.append(f'{dep_var}_residual_rank')
            phenotypes = data[:, col_names.index(f'{dep_var}_residual')]
        else:
            col_names.append(f'{dep_var}_rank')
            phenotypes = data[:, col_names.index(dep_var)]
        not_nan = ~np.isnan(phenotypes)
        phenotypes = phenotypes[not_nan]
        sort = np.argsort(phenotypes)
        ranks = np.empty_like(sort)
        ranks[sort] = np.arange(phenotypes.shape[0])
        data[not_nan, -1] = ranks
    return data


def inverse_normalize_ranks(
        readme,
        data,
        col_names,
        dep_vars,
        residuals=True):
    """
    col_names will be modified
    """
    readme.write("Inverse normalizing phenotype residuals ranks to the standard normal distribution "
                 "via the transformation rank -> normal_quantile((rank + 0.5)/nsamples).\n")
    for dep_var in dep_vars:
        data = np.concatenate((
            data,
            np.full((data.shape[0], 1), np.nan)
        ), axis=1)
        if residuals:
            col_names.append(f'{dep_var}_inv_norm_rank_residual')
            ranks = data[:, col_names.index(f'{dep_var}_residual_rank')]
        else:
            col_names.append(f'{dep_var}_inv_norm_rank')
            ranks = data[:, col_names.index(f'{dep_var}_rank')]
        not_nan = ~np.isnan(ranks)
        ranks = ranks[not_nan]
        data[not_nan, -1] = scipy.stats.norm.ppf((ranks + 0.5)/ranks.shape[0])
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
        '--residual-model',
        help='which model to use for calculating residuals when regressing out covaraites',
        type=str,
        choices=['linear', 'krr', 'random-forest', 'adaboost'],
        default='linear'
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
        readme.write(f"Subsetting to samples from file {floc}. ")
        sample_subset = np.genfromtxt(
            floc,
            usecols=[0],
            skip_header=1,
            delimiter=" "
        ).reshape(-1, 1)
        data = merge_arrays(sample_subset, data)
        readme.write(f"{data.shape[0]} remaining samples\n")
        readme.flush()

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
            if args.residual_model == 'linear':
                data = get_residuals_linear(
                    readme,
                    assoc_dir,
                    data,
                    col_names,
                    dep_vars,
                    indep_vars,
                    pheno_indep_vars
                )
            elif args.residual_model == 'krr':
                data = get_residuals_kernel_ridge(
                    readme,
                    assoc_dir,
                    data,
                    col_names,
                    dep_vars,
                    indep_vars,
                    pheno_indep_vars
                )
            elif args.residual_model == 'random-forest':
                data = get_residuals_random_forest(
                    readme,
                    assoc_dir,
                    data,
                    col_names,
                    dep_vars,
                    indep_vars,
                    pheno_indep_vars
                )
            elif args.residual_model == 'adaboost':
                data = get_residuals_adaboost(
                    readme,
                    assoc_dir,
                    data,
                    col_names,
                    dep_vars,
                    indep_vars,
                    pheno_indep_vars
                )

        data = rank_phenotypes(readme, data, col_names, dep_vars, residuals=not args.plink)
        data = inverse_normalize_ranks(readme, data, col_names, dep_vars, residuals=not args.plink)

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
                "Performing linear association of inverse normal transformed ranks of residual phenotypes "
                "against 0 to 2 genotype dosage of alt (== non-ref) allele\n"
            )
            readme.write(
                "Not filtering loci for INFO or calls for prob.\n"
                "Filtering loci roughly equivalent to minor allele count < 20. "
                "Specifically, treating 0.5 < alt dosage 1.5 as 1 alt allele and "
                "alt dosage > 1.5 as 2 alt alleles and calculating minor allele "
                "count that way\n."
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
                "Running linear associations of dosage vs inverse normal transformed phenotypes "
                "with plink.\n"
                "Handing to plink all relevant covariates.\n"
                "Filtering loci with <20 nonmajor alleles as recommended by plink docs: "
                "https://www.cog-genomics.org/plink/2.0/assoc\n"
            )
            readme.flush()
            indep_vars.remove(col_names.index('intercept'))
            run_plink_associations(assoc_dir, dep_vars, indep_vars, pheno_indep_vars, col_names,
                                   chromosome=args.chromosome)
        print("Association run complete.")

if __name__ == "__main__":
    main()

