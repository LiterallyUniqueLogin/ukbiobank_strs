#!/usr/bin/env python3

import argparse
import datetime
import sys

import numpy as np

import python_array_utils as utils

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
parser.add_argument('data_fname')
parser.add_argument('outprefix')
parser.add_argument('samples')
parser.add_argument('assessment_ages_npy')
parser.add_argument(
    '--categorical-covar-names',
    nargs='*',
    default=[]
)
parser.add_argument(
    '--categorical-covar-files',
    nargs='*',
    default=[]
)

args = parser.parse_args()

outprefix = args.outprefix

assert len(args.categorical_covar_names) == len(args.categorical_covar_files)

with open(f'{outprefix}_README.txt', 'w') as readme, \
        open(f'{outprefix}_covar_names.txt', 'w') as covar_names:
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    data_fname = args.data_fname
    readme.write(f"Run date: {today}\n")
    readme.write(
        f"Loading from txt file "
        f" {data_fname} \n"
    )

    data = np.genfromtxt(
        data_fname,
        skip_header = 1,
        delimiter='\t'
    )[:, 1:-1]
    # drop first and last rows which because of the way data is extracted
    # and then read by numpy are always nans

    readme.write(f"Subsetting to samples at {args.samples}\n")
    # load samples that have passed qc
    samples = np.genfromtxt(
        args.samples,
        skip_header=1
    )

    # subset to those samples
    data = utils.merge_arrays(
        samples.reshape(-1, 1),
        data
    )

    # number of sample with this phenotype at any assessment
    num_samples = np.sum(np.any(~np.isnan(data[:, 1:]), axis=1))

    # drop samples with categorical covars with less than this
    # number of samples or fraction of samples
    cat_drop_num = 50
    cat_drop_frac = 0.001 #0.1%

    covar_datas = []
    reverse_covar_hashes = []
    for covar_name, covar_fname in zip(args.categorical_covar_names, args.categorical_covar_files):
        readme.write(
            f"Loading categorical covar with name {covar_name} from txt file "
            f" {covar_fname}\n"
        )
        obj_covar_data = np.genfromtxt(
            covar_fname,
            skip_header = 1,
            delimiter='\t',
            dtype=object
        )[:, 1:-1]
        assert len(obj_covar_data.shape) == len(data.shape)

        # do some work to hash the array so it's faster to merge
        # (object arrays are slow to work with)
        hash_covar_data = np.full(obj_covar_data.shape, np.nan, dtype=float)
        hash_covar_data[:, 0] = obj_covar_data[:, 0].astype(int) #participant ids
        obj_covar_data = obj_covar_data[:, 1:]
        unique_vals = np.unique(obj_covar_data[obj_covar_data != b''])
        unique_hashes = np.unique([hash(val) for val in unique_vals])
        assert len(unique_hashes) == len(unique_vals)
        hash_view = hash_covar_data[:, 1:]
        reverse_covar_hashes.append({})
        for val in unique_vals:
            hash_view[obj_covar_data == val] = float(hash(val))
            reverse_covar_hashes[-1][float(hash(val))] = val.decode()
        assert len(
            np.unique(hash_covar_data[:, 1:][~np.isnan(hash_covar_data[:, 1:])])
        ) == len(unique_vals)

        covar_datas.append(utils.merge_arrays(
            data[:, 0:1], #sample_ids in order
            hash_covar_data
        ))
        assert covar_datas[-1].shape == data.shape

    ages = np.load(args.assessment_ages_npy)
    ages = utils.merge_arrays(data[:, 0:1], ages)
    covar_names.write('age\n')
    readme.write(
        "Choosing phenotype value and age for each participant based on the first "
        "visit for which this phenotype had a recorded value. If the participant "
        "had this phenotype measured at multiple visits, only the phenotype value and age "
        "at the first visit are being used. The age "
        "is being loaded from the shared_covars file "
        f"{args.assessment_ages_npy} . An additional "
        "dummy covariate indicating visit number is being added "
        "for each visit whose data is used beyond the first.\n"
    )

    # figure out what assessments the data was recorded at
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

    # move data into visit_aligned_data as appropriate
    # 3 columns: ID, phenotype, age
    # one additional column to be added as a dummy variable for each
    # assessment beyond the first
    visit_aligned_data = np.full((data.shape[0], 3), np.nan)
    visit_aligned_data[:, 0] = data[:, 0]

    # move each covar into visit_aligned_covars likewise
    # at the end those will be turned into categorical dummy variables
    # these are just 1D arrays with the same ordering as visit_aligned_data
    visit_aligned_covars = [
        np.full((data.shape[0], 1), np.nan) for _ in covar_datas
    ]

    # set phenotypes, then ages and other covars
    visit_aligned_data[:, 1] = data[:, assess_nums[0] + 1]
    visit_aligned_data[:, 2] = ages[:, assess_nums[0] + 1]

    for (covar_data, visit_aligned_covar) in zip(covar_datas, visit_aligned_covars):
        visit_aligned_covar[:, 0] = covar_data[:, assess_nums[0] + 1]

    # assert that the first recording of data shouldn't be dropped
    num_first_rec = np.sum(~np.isnan(data[:, assess_nums[0] + 1]))
    assert num_first_rec >= cat_drop_num
    assert num_first_rec >= num_samples*cat_drop_frac

    # now add samples from additional assessments
    for assess_num in assess_nums[1:]:
        nans = np.isnan(visit_aligned_data[:, 1])
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
        visit_aligned_data[nans, 1] = data[nans, assess_num + 1]
        visit_aligned_data[nans, 2] = ages[nans, assess_num + 1]
        for covar_data, visit_aligned_covar in zip(covar_datas, visit_aligned_covars):
            visit_aligned_covar[nans, 0] = covar_data[nans, assess_num + 1]

        filled_lines = nans & ~np.isnan(visit_aligned_data[:, 1])
        # add a dummy coavriate for this visit as a batch effect
        if np.any(filled_lines):
            covar_names.write(
                'pheno_from_' + reverse_assessment_dict[assess_num].replace('-','_') + "\n"
            )
            visit_aligned_data = np.concatenate(
                (visit_aligned_data, filled_lines.reshape(-1, 1)),
                axis=1
            )

    data = visit_aligned_data
    # we will subset to those samples which only have measured data
    measured_samples = ~np.any(np.isnan(data), axis=1)
    data = data[measured_samples, :]

    for covar_name, covar_data, reverse_covar_hash in \
            zip(args.categorical_covar_names, visit_aligned_covars, reverse_covar_hashes):
        unmeasured_covar_data = covar_data[~measured_samples, :]
        covar_data = covar_data[measured_samples, :]
        assert covar_data.shape[1] == 1
        categories, counts = np.unique(covar_data[~np.isnan(covar_data)], return_counts=True)
        cat_names = [reverse_covar_hash[cat] for cat in categories]
        max_cat_idx = np.argmax(counts)
        assert len(categories) <= 10
        using_covar = False
        for category, cat_name in reverse_covar_hash.items():
            if cat_name not in cat_names:
                if not np.any(unmeasured_covar_data == category):
                    readme.write(
                        f"No qc'ed samples of this ethnicity had value {cat_name} "
                        f'for categorical covariate {covar_name}.\n'
                    )
                else:
                    readme.write(
                        f"No qc'ed samples of this ethnicity for which this phenotype "
                        f'was measured had value {cat_name} '
                        f'for categorical covariate {covar_name}.\n'
                    )
                continue

            if cat_name == cat_names[max_cat_idx]:
                continue
            cat_idx = cat_names.index(cat_name)
            cat_count = counts[cat_idx]
            if cat_count <= cat_drop_num:
                readme.write(
                    f"Dropping {cat_count} samples with value {cat_name} for categorical "
                    f"covariate {covar_name} because that is less than "
                    f"{cat_drop_num} total samples.\n"
                )
                continue
            elif cat_count <= num_samples*cat_drop_frac:
                readme.write(
                    f"Dropping {cat_count} samples with value {cat_name} for categorical "
                    f"covariate {covar_name} because that is less than the fraction "
                    f"{cat_drop_frac} of the total samples for this phenotype "
                    f"which is {num_samples}.\n"
                )
                continue
            using_covar = True
            covar_names.write(f"{covar_name}_is_{cat_name}\n")
            data = np.concatenate((data, covar_data == category), axis=1)
        if using_covar:
            readme.write(
                f'Using categorical covar {covar_name}, default category when all dummy '
                f'variables are False is {cat_names[max_cat_idx]}\n'
            )
        else:
            readme.write(
                f'Excluding categorical covar {covar_name} '
                '- all categories after the first were '
                'already excluded due to small sizes\n'
            )

    # drop samples which are missing the phenotype
    data = data[~np.isnan(data[:, 1]), :]
    # drop samples which are missing a covariate
    has_all_covars = np.all(~np.isnan(data), axis=1)
    if not np.all(has_all_covars):
        readme.write(
            f'Dropping {np.sum(~has_all_covars)} samples that are missing '
            'covariates. Data:\n'
        )
        with np.printoptions(threshold=sys.maxsize):
            readme.write(str(data[~has_all_covars, :]))
        data = data[has_all_covars, :]
    else:
        readme.write('All samples have all covariates.\n')

    np.save(f'{outprefix}.npy', data)

