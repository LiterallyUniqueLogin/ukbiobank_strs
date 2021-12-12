#!/usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--ethnicity', default='white_brits')

ethnicity = parser.parse_args().ethnicity

ukb = os.environ['UKB']



with open(f'{ukb}/microarray/ukb46122_hap_chr1_v2_s487314.sample') as sample_file:
    next(sample_file) # header
    next(sample_file) # sample ID 0
    samples = set(line.split()[0] for line in sample_file if line[0] != '-')

print(f'# samples with genetic data that have not withdrawn before microarray v2 publication: {len(samples)}')

with open(f'{ukb}/sample_qc/common_filters/remove/withdrawn.sample') as withdrawn_file:
    next(withdrawn_file) # header
    samples.difference_update(set(line.strip() for line in withdrawn_file))

not_withdrawn_total = len(samples)
print(f'# samples with genetic data that have not withdrawn up to now: {not_withdrawn_total}')

with open(f'{ukb}/sample_qc/common_filters/ethnicity/{ethnicity}.sample') as ethnicity_file:
    next(ethnicity_file) # header
    samples.intersection_update(set(line.split()[0] for line in ethnicity_file))

ethnic_subset_total = len(samples)
print(f'# {ethnicity} samples with genetic data that have not withdrawn up to now: {ethnic_subset_total} (%{ethnic_subset_total/not_withdrawn_total * 100})')

with open(f'{ukb}/sample_qc/common_filters/remove/sex_mismatch.sample') as sex_mismatch_file:
    next(sex_mismatch_file) # header
    sex_mismatch = samples.intersection(set(line.split()[0] for line in sex_mismatch_file))

print(f'# sex mismatch {ethnicity} samples with genetic data that have not withdrawn up to now: {len(sex_mismatch)}')

with open(f'{ukb}/sample_qc/common_filters/remove/sex_aneuploidy.sample') as sex_aneuploidy_file:
    next(sex_aneuploidy_file) # header
    sex_aneuploidy = samples.intersection(set(line.split()[0] for line in sex_aneuploidy_file))
    sex_aneuploidy.difference_update(sex_mismatch)

print(f'# additional sex aneuploidy {ethnicity} samples with genetic data that have not withdrawn up to now: {len(sex_aneuploidy)}')

with open(f'{ukb}/sample_qc/common_filters/remove/low_quality.sample') as low_quality_file:
    next(low_quality_file) # header
    low_quality = samples.intersection(set(line.split()[0] for line in low_quality_file))
    low_quality_and_sex_aneuploidy = low_quality.intersection(sex_aneuploidy)
    low_quality_and_sex_mismatch = low_quality.intersection(sex_mismatch)

print(f'# low quality {ethnicity} samples with genetic data that have not withdrawn up to now: {len(low_quality)}')
print(f'# low quality and sex aneuploidy {ethnicity} samples with genetic data that have not withdrawn up to now: {len(low_quality_and_sex_aneuploidy)}')
print(f'# low quality and sex mismatch {ethnicity} samples with genetic data that have not withdrawn up to now: {len(low_quality_and_sex_mismatch)}')

samples.difference_update(low_quality)
samples.difference_update(sex_aneuploidy)
samples.difference_update(sex_mismatch)

print(f'# QCed {ethnicity} samples with genetic data that have not withdrawn up to now: {len(samples)} ({len(samples)/ethnic_subset_total*100}% of {ethnicity} subset total)')
