#!/usr/bin/env python3

import argparse
import os
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
import sklearn.model_selection
import sklearn.neighbors

import python_array_utils as utils

ukb = os.environ['UKB']

def fit_kde(data: np.ndarray,
            random_state: int = 13,
            max_train = int(1e4)) -> sklearn.neighbors.KernelDensity:
    """
    Parameters
    ----------
    data:
        2D array of shape nsamples x nfeatures
    """
    # only train on up to 10k loci for speed. Select a random subset.
    # it would be nice to train on the same subset of loci for each strata
    # but that isn't feasible if some strata have many more nan's than
    # others
    assert data.shape[0] > 1000
    assert not np.any(np.isnan(data))

    # code from
    # https://jakevdp.github.io/PythonDataScienceHandbook/05.13-kernel-density-estimation.html

    # Use gridsearch to choose the bandwidth
    kfold = sklearn.model_selection.ShuffleSplit(
        n_splits=5,
        train_size=min(max_train, int(np.floor(data.shape[0]*.85))),
        random_state=random_state
    )
    bandwidths = 10 ** np.linspace(-1.8, 1.8, 40)
    grid = sklearn.model_selection.GridSearchCV(
        sklearn.neighbors.KernelDensity(kernel='gaussian'),
        {'bandwidth': bandwidths},
        cv=kfold,
        n_jobs=5
    )
    grid.fit(data)
    bandwidth = grid.best_params_['bandwidth']

    #compute the kde with the best bandwidth
    kde = sklearn.neighbors.KernelDensity(kernel='gaussian',
                                          bandwidth=bandwidth)
    kde.fit(data)
    return kde


# copied from plot_stats branch of trtools
def plot_1D_kde(
        data: np.ndarray,
        xlabel: str,
        title: str,
        fname: str,
        strata_labels: Optional[Dict[float, str]] = None):
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
        or a 2D array with two columns where the first column
        should be stratified by the second
    xlabel:
        the x label for the graph
    title:
        the title for the graph
    fname:
        the file name to save the graph. Must include the extension,
        and one that matplotlib will recognize so that it produces
        a file of that type.
    strata_labels:
        if data is 2D, then a dict from strata values to labels
    """
    if len(data.shape) == 1:
        data = data.reshape(-1, 1)
    else:
        assert len(data.shape) == 2
        assert data.shape[1] == 2
        assert strata_labels is not None
        strata_values = np.unique(data[:, 1])
        assert len(strata_labels) == len(strata_values)
        for stratum_value in strata_values:
            assert stratum_value in strata_labels

    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Probability density")

    if data.shape[1] == 1:
        strata = {None: data}
    else:
        strata = {}
        for stratum_value in strata_values:
            stratum = data[data[:, 1] == stratum_value, 0]
            stratum = stratum[~np.isnan(stratum)]
            stratum = stratum.reshape(-1, 1)
            strata[stratum_value] = stratum

    # Fit and plot each stratum individually
    for stratum_value, stratum in strata.items():
        kde = fit_kde(stratum)

        min_val = np.min(stratum)
        max_val = np.max(stratum)
        eps = (max_val - min_val)/10e3
        xs = np.arange(min_val - eps, max_val + eps, eps)
        curve = np.exp(kde.score_samples(xs.reshape(-1, 1)))

        # plot
        ax.fill_between(xs, curve, alpha=0.5)

    if data.shape[1] > 1:
        ax.legend()
    plt.savefig(fname)


def plot_2D_kde(
        data: np.ndarray,
        xlabel: str,
        ylabel: str,
        title: str,
        fname: str):
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
        A 2D array with two columns where the first column is the
        x variable and the second column is the y variable.
    xlabel:
        the x label for the graph
    ylabel:
        the y label for the graph
    title:
        the title for the graph
    fname:
        the file name to save the graph. Must include the extension,
        and one that matplotlib will recognize so that it produces
        a file of that type.
    """
    assert len(data.shape) == 2
    assert data.shape[1] == 2
    assert not np.any(np.isnan(data))

    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    kde = fit_kde(data, max_train = 50)

    min_xval = np.min(data[:, 0])
    max_xval = np.max(data[:, 0])
    xeps = (max_xval - min_xval)/10e3
    xs_1d = np.arange(min_xval - xeps, max_xval + xeps, xeps)

    min_yval = np.min(data[:, 0])
    max_yval = np.max(data[:, 0])
    yeps = (max_yval - min_yval)/10e3
    ys_1d = np.arange(min_yval - yeps, max_yval + yeps, yeps)

    xs, ys = np.meshgrid(xs_1d, ys_1d)
    xy = np.vstack([ys.ravel(), xs.ravel()]).T

    probs = np.exp(kde.score_samples(xy))
    probs = probs.reshape(xy.shape)

    cs = ax.contourf(xs, ys, probs, colors=plt.cm.Blues)
    cb = fig.colorbar(cs)
    cb.ax.set_title("Probability density")

    plt.savefig(fname)


def plot_histogram(
        data: np.ndarray,
        xlabel: str,
        title: str,
        fname: str,
        strata_labels: Dict[float, str]):
    assert not np.any(np.isnan(data))

    min_val = np.min(data[:, 0])
    max_val = np.max(data[:, 0])
    unique_points = np.unique(data[:, 0])
    print(unique_points)
    nbins = len(unique_points)
    while True:
        if len(np.unique(np.floor(
            (unique_points - min_val)*nbins/(max_val - min_val)
        ))) == len(unique_points):
            break
        nbins += 1

    fig, ax = plt.subplots()

    for stratum_value, stratum_label in strata_labels.items():
        ax.hist(
            data[data[:, 1] == stratum_value, 0],
            nbins,
            range = (min_val, max_val),
            density = True,
            label = stratum_label,
            alpha=0.5
        )

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Probability density")
    ax.legend()

    plt.savefig(fname)


def plot_phenotype_by_sex(phenotype):
    shared_covars = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')
    pheno_data = np.load(f'{ukb}/traits/phenotypes/{phenotype}.npy')
    with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt') as unit_file:
        unit = next(unit_file).strip()

    data = utils.merge_arrays(shared_covars, pheno_data)
    data = data[:, [shared_covars.shape[1], 1]]
    data = data[~np.any(np.isnan(data), axis=1), :]

    if len(np.unique(data[:, 0])) < 2000:
        plot_histogram(
            data,
            f'{phenotype} ({unit})',
            phenotype.capitalize() + ' x Sex distribution',
            f'{ukb}/traits/phenotypes/{phenotype}_distribution_by_sex.png',
            {1: 'male', 2: 'female'}
        )
    else:
        plot_1D_kde(
            data,
            f'{phenotype} ({unit})',
            phenotype.capitalize() + ' x Sex distribution',
            f'{ukb}/traits/phenotypes/{phenotype}_distribution_by_sex.png',
            {1: 'male', 2: 'female'}
        )


def plot_phenotype_by_age(phenotype):
    pheno_data = np.load(f'{ukb}/traits/phenotypes/{phenotype}.npy')
    with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt') as unit_file:
        unit = next(unit_file).strip()
    pheno_data = pheno_data[:, [1, 2]]
    pheno_data = pheno_data[~np.any(np.isnan(pheno_data), axis=1), :]
    plot_2D_kde(
        pheno_data,
        f'{phenotype} ({unit})',
        'age (years)',
        phenotype.capitalize() + ' x Age distribution',
        f'{ukb}/traits/phenotypes/{phenotype}_distribution_by_age.png'
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype', type=str)
    parser.add_argument('covar', choices=['sex', 'age'])
    args = parser.parse_args()

    if args.covar == 'sex':
        plot_phenotype_by_sex(args.phenotype)
    else:
        plot_phenotype_by_age(args.phenotype)

