#!/usr/bin/env python3

import argparse
import datetime
import os

import numpy as np

import python_array_utils as utils

ukb = os.environ['UKB']

assessment_dict = {
    'init-assessment': 0,
    'repeat-assessment-1': 1,
    'imaging-visit': 2,
    'repeat-imaging-1': 3
}
reverse_assessment_dict = {}
for key, value in assessment_dict.items():
    reverse_assessment_dict[value] = key

parser = argparse.ArgumentParser()
parser.add_argument('phenotype_name')
parser.add_argument('phenotype_field_id')
parser.add_argument('unit')
parser.add_argument(
    'age',
    choices = set(assessment_dict.keys()).union({'first-available'})
)
parser.add_argument(
    '--instance-id'
)
parser.add_argument(
    '--categorical-covars',
    nargs='+',
    default=[]
)

args = parser.parse_args()

assert args.unit != 'binary'

phenotype = args.phenotype_name

with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt', 'w') as unit_file:
    unit_file.write(f'{args.unit}\n')

with open(f'{ukb}/traits/phenotypes/{phenotype}_README.txt', 'w') as readme, \
        open(f'{ukb}/traits/phenotypes/{phenotype}_covar_names.txt', 'w') as covar_names:
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    data_fname = f'{ukb}/main_dataset/extracted_data/{phenotype}_{args.phenotype_field_id}.txt'
    readme.write(f"Run date: {today}\n")
    readme.write(
        f"Loading phenotype {phenotype} from txt file "
        f" {data_fname} \n"
    )

    data = np.genfromtxt(
        data_fname,
        skip_header = 1,
        delimiter='\t'
    )[:, 1:-1]
    # drop first and last rows which because of the way data is extracted
    # and then read by numpy are always nans

    # number of sample with this phenotype at any assessment
    num_samples = np.sum(np.any(~np.isnan(data[:, 1:]), axis=1))

    # drop samples with categorical covars with less than this
    # number of samples or fraction of samples
    cat_drop_num = 50
    cat_drop_frac = 0.001 #0.1%

    covar_datas = []
    for covar in args.categorical_covars:
        covar_name, covar_id = covar.split(',')
        covar_fname = f'{ukb}/main_dataset/extracted_data/{covar_name}_{covar_id}.txt'
        readme.write(
            f"Loading categorical covar with field id {covar_id} from txt file "
            f" {covar_fname}\n"
        )
        covar_datas.append(utils.merge_arrays(
            data[:, 0:1], #sample_ids in order
            np.genfromtxt(
                covar_fname,
                skip_header = 1,
                delimiter='\t'
            )[:, 1:-1])
        )
        assert covar_datas[-1].shape == data.shape

    ages = np.load(f'{ukb}/traits/shared_covars/assessment_ages.npy')
    if args.age in assessment_dict:
        covar_names.write('age\n')
        readme.write(
            f"Choosing age for each participant corresponding to the visit "
            f"'{args.age}'. This is being loaded from the shared_covars file "
            f"{ukb}/traits/shared_covars/assessment_ages.npy\n"
        )
        # assuming data was only taken at that assessment, not that it was
        # taken at multiple and we only want one
        assert data.shape[1] == 2
        assess_aligned_covar_datas = covar_datas

        col = assessment_dict[args.age]

        data = utils.merge_arrays(data, ages[:, [0, col + 1]])

    elif args.age == 'first-available':
        covar_names.write('age\n')
        readme.write(
            "Choosing phenotype value and age for each participant based on the first "
            "visit for which this phenotype had a recorded value. If the participant "
            "had this phenotype measured at multiple visits, only the phenotype value and age "
            "at the first visit are being used. The age "
            "is being loaded from the shared_covars file "
            f"{ukb}/traits/shared_covars/asssessment_ages.npy . An additional "
            "dummy covariate indicating visit number is being added "
            "for each visit whose data is used beyond the first.\n"
        )
        with open(data_fname) as data_file:
            data_header = next(data_file)
        col_names = data_header.split()
        assert col_names[0] == 'eid'
        field_id = col_names[1].split('-')[0]
        assess_nums = []
        for idx, col_name in enumerate(col_names[1:]):
            assert len(col_name) == len(field_id) + 4
            assert col_name.startswith(field_id)
            assert col_name.endswith('.0')
            assess_num = int(col_name[-3])
            assess_nums.append(assess_num)

        assess_aligned_covar_datas = []
        for _ in covar_datas:
            assess_aligned_covar_data = np.full((data.shape[0], 2), np.nan)
            assess_aligned_covar_data[:, 0] = data[:, 0]
            assess_aligned_covar_datas.append(assess_aligned_covar_data)

        ages = utils.merge_arrays(data[:, 0:1], ages)
        # move data into new_data as appropriate
        # 3 initial columns: ID, phenotype, age
        # one additional column to be added as a dummy variable for each
        # assessment beyond the first
        new_data = np.full((data.shape[0], 3), np.nan)
        new_data[:, 0] = data[:, 0]
        # set phenotypes, then ages
        new_data[:, 1] = data[:, assess_nums[0] + 1]
        new_data[:, 2] = ages[:, assess_nums[0] + 1]
        for covar_data, assess_aligned_covar_data in zip(covar_datas, assess_aligned_covar_datas):
            assess_aligned_covar_data[:, 1] = covar_data[:, assess_nums[0] + 1]
        # assert that the first recording of data shouldn't be dropped
        num_first_rec = np.sum(~np.isnan(data[:, assess_nums[0] + 1]))
        assert num_first_rec >= cat_drop_num
        assert num_first_rec >= num_samples*cat_drop_frac
        for assess_num in assess_nums[1:]:
            nans = np.isnan(new_data[:, 1])
            assessment_data = data[nans, assess_num + 1]
            num_assessment_data = np.sum(~np.isnan(assessment_data))
            if num_assessment_data <= cat_drop_num:
                readme.write(
                    f"Dropping {num_assessment_data} samples whose first "
                    f"(non yet dropped) data recording for this phenotype occurred on assessment "
                    f"{assess_num} because that is less than "
                    f"{cat_drop_num} total samples.\n"
                )
                continue
            elif num_assessment_data <= num_samples*cat_drop_frac:
                readme.write(
                    f"Dropping {num_assessment_data} samples whose first "
                    f"(non yet dropped) data recording for this phenotype occurred on assessment "
                    f"{assess_num} because that is less than the fraction "
                    f"{cat_drop_frac} of the total samples for this phenotype "
                    f"which is {num_samples}.\n"
                )
                continue
            new_data[nans, 1] = data[nans, assess_num + 1]
            new_data[nans, 2] = ages[nans, assess_num + 1]
            for covar_data, assess_aligned_covar_data in zip(covar_datas, assess_aligned_covar_datas):
                assess_aligned_covar_data[nans, 1] = covar_data[nans, assess_nums[0] + 1]
            filled_lines = nans & ~np.isnan(new_data[:, 1])
            if np.any(filled_lines):
                covar_names.write(
                    'pheno_from_' + reverse_assessment_dict[assess_num].replace('-','_') + "\n"
                )
                new_data = np.concatenate((new_data, filled_lines.reshape(-1, 1)), axis=1)
        data = new_data
    else:
        raise Exception("Age option not understood")
        # currently don't have a setting for age is None
        # because wouldn't know how to handle categorical covars
        # or how to merge across different timestamps

    for covar, covar_data in zip(args.categorical_covars, assess_aligned_covar_datas):
        assert covar_data.shape[1] == 2
        covar_name = covar.split(',')
        cats = np.unique(covar_data[~np.isnan(covar_data[:, 1]), 1])
        assert len(cats) <= 10
        n_cat_0  = np.sum(covar_data == cats[0])
        assert n_cat_0 > cat_drop_num
        assert n_cat_0 > num_samples*cat_drop_frac
        for cat in cats[1:]:
            is_cat = covar_data[:, 1] == cat
            n_cat = np.sum(is_cat)
            if n_cat <= cat_drop_num:
                readme.write(
                    f"Dropping {n_cat} samples with value {cat} for categorical "
                    f"covariate {covar} because that is less than "
                    f"{cat_drop_num} total samples.\n"
                )
                continue
            elif n_cat <= num_samples*cat_drop_frac:
                readme.write(
                    f"Dropping {n_cat} samples with value {cat} for categorical "
                    f"covariate {covar} because that is less than the fraction "
                    f"{cat_drop_frac} of the total samples for this phenotype "
                    f"which is {num_samples}.\n"
                )
                continue
            covar_names.write("{covar_name}_is_{cat}\n")
            data = np.concatenate((data, is_cat), axis=1)

    # drop samples which are missing the phenotype
    data = data[~np.isnan(data[:, 1]), :]
    # drop samples which are missing a covariate
    has_all_covars = np.all(~np.isnan(data), axis=1)
    readme.write(
        f'Dropping {np.sum(~has_all_covars)} samples that are missing '
        'covariates\n.'
    )
    data = data[has_all_covars, :]

    np.save(f'{ukb}/traits/phenotypes/{phenotype}.npy', data)

