#!/usr/bin/env python3

import argparse
import json
import os

ukb = os.environ['UKB']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('readme_fnames_json_fname')
    parser.add_argument('samples_fnames_json_fname')
    args = parser.parse_args()
    with open(args.readme_fnames_json_fname) as f:
        readme_fnames_dict = json.loads(next(f))
    with open(args.samples_fnames_json_fname) as f:
        sample_fnames_dict = json.loads(next(f))

    assert set(readme_fnames_dict) == set(sample_fnames_dict)
    pairs = sorted(readme_fnames_dict)
    with open(f'{ukb}/export_scripts/results/ethnicity_associations.tab', 'w') as outfile:
        outfile.write(
            'ethnicity\t'
            'phenotype\t'
            'n_individuals_tested\t'
            'measurement_visits\t'
            'categorical_covariates\n'
        )
        for pair in pairs:
            ethnicity, phenotype = pair.split('__')
            outfile.write(f'{ethnicity}\t{phenotype}\t')
            with open(sample_fnames_dict[pair]) as f:
                n_samples = len(f.readlines()) - 1
                outfile.write(f'{n_samples}\t')
            with open(readme_fnames_dict[pair]) as readme:
                text = readme.read()
                if 'assessment 1' in text:
                    outfile.write('initial_assessment\t')
                else:
                    outfile.write('initial_assessment,first_repeat_assessment\t')

                if 'device_id' in text:
                    outfile.write('device_id\n')
                elif 'aliquot' in text:
                    values = []
                    for value in range(0, 5):
                        if f'value {value}' not in text:
                            values.append(value)
                    values = tuple(values)
                    if len(values) > 1:
                        outfile.write(f'aliquot, {values}\n')
                    else:
                        outfile.write('None\n')
                else:
                    outfile.write('None\n')

    with open(f'{ukb}/export_scripts/results/ethnicity_associations_README.txt', 'w') as readme:
        readme.write('ethnicity: the ethnic group results were replicated in\n')
        readme.write('phenotype: the name of the phenotype\n')
        readme.write(
            'n_individuals_tested: the number of individuals included in the association testing '
            'i.e. unrelated high-quality samples with the phenotype\n'
        )
        readme.write(
            'measurement_visits: the different visits measurements were drawn '
            '(always choosing the first if measurements were taken at multiple visits). '
            'If more than one, then a categorical covariate is added for each beyond the first. '
        )
        readme.write(
            'categorical_covariates: categorical covariates included in the '
            'association tests other than those for visits (see above). Covariates that all aassume all values '
            'are listed by name, covariates that assume only some values are listed as pairs (covariate name, '
            'tuple of assumed values)\n'
        )

if __name__ == '__main__':
    main()

