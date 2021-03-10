#!/usr/bin/env python3

import argparse
import datetime
import os
import os.path

import numpy as np
import sklearn.model_selection
import statsmodels.regression.linear_model

import python_array_utils as utils

ukb = os.environ['UKB']

def get_linear_residuals(readme, data, covar_names):
    splitter = sklearn.model_selection.ShuffleSplit(
        n_splits = 5,
        train_size = 0.9,
        random_state = hash('linear') % 2**32
    )
    n_samples = data.shape[0]
    readme.write(
        f"Using 5 splits with 90% of the data ({0.9*n_samples:.0f} samples) "
        f"as training and 10% of the data ({0.1*n_samples:.0f} samples) as validation "
        "in order to calculate average validation RMSE for the regression.\n"
    )
    readme.flush()
    y = data[:, 1]
    temp_covars = [covar_names.index('age') + 2, covar_names.index('sex') + 2]
    X = data[:, temp_covars]
    #X = data[:, 2:]
    total_rmse = 0
    for train, test in splitter.split(X):
        model = statsmodels.regression.linear_model.OLS(y[train], X[train, :])
        reg = model.fit()
        total_rmse += np.sqrt(np.mean((y[test] - reg.predict(X[test, :]))**2))
    avg_rmse = total_rmse/splitter.get_n_splits()
    readme.write(f"Average RMSE: {avg_rmse:.4f}\n")
    readme.flush()

    model = statsmodels.regression.linear_model.OLS(y, X)
    reg = model.fit()
    readme.write("Per covaraite results:\n")
    readme.write("indep_var\tp\tcoeff\n")
    #for idx, var in enumerate(covar_names):
    for idx, var in enumerate(['age', 'sex']):
        readme.write(f'{var}\t{reg.pvalues[idx]:.2e}\t'
                     f'{reg.params[idx]}\n')
    readme.flush()
    return y - reg.predict(X)

def main():  # noqa: D103
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'regression_model',
        type=str,
        choices=['linear']
    )
    parser.add_argument('phenotype')
    args = parser.parse_args()

    model = args.regression_model
    phenotype = args.phenotype

    with open(f'{ukb}/traits/adjusted_srin_phenotypes/{phenotype}_{model}_README.txt', 'w') as readme:
        today = datetime.datetime.now().strftime("%Y_%m_%d")
        readme.write(f"Run date: {today}\n")
        readme.flush()

        data = np.load(f'{ukb}/traits/subset_rin_phenotypes/{phenotype}.npy')
        with open(f'{ukb}/traits/phenotypes/{phenotype}_covar_names.txt') as pheno_covar_file:
            covar_names = [line.strip() for line in pheno_covar_file if line.strip() != '']
        # the last three columns are assessment dates which aren't in and of themselves covariates
        covars = np.load(f'{ukb}/traits/shared_covars/shared_covars.npy')[:, :-3]
        with open(f'{ukb}/traits/shared_covars/covar_names.txt') as covar_file:
            covar_names.extend([line.strip() for line in covar_file if line.strip() != ''])

        data = utils.merge_arrays(data, covars)
        assert not np.any(np.isnan(data))
        # now first col is IDs, second col is phenotype, all other cols are covars
        assert data.shape[1] == len(covar_names) + 2
        # data should be standard normally distributed
        assert np.isclose(np.mean(data[:, 1]), 0)
        assert np.isclose(np.std(data[:, 1]), 1)

        if model == 'linear':
            residuals = get_linear_residuals(readme, data, covar_names)
        else:
            raise ValueError('Unimplemented model type')

        assert not np.any(np.isnan(residuals))
        np.save(
            f'{ukb}/traits/adjusted_srin_phenotypes/{phenotype}_{model}.npy',
            np.concatenate((data[:, 0:1], residuals.reshape(-1, 1)), axis=1)
        )

if __name__ == "__main__":
    main()

