#!/usr/bin/env python3
import argparse
import os
import sys
from typing import Dict, Set

parser = argparse.ArgumentParser()
parser.add_argument(
    'run_name',
    help=('Verifies the unrelatedness of the sample file '
          '$UKB/sample_qc/runs/{run_name}/combined_unrelated.sample')
)
args = parser.parse_args()
run_name = args.run_name
ukb = os.environ['UKB']

def error(msg):
    print(msg, file=sys.stderr)
    sys.exit(-1)

def load_samples(fname):
    samples = set()
    with open(fname) as sample_file:
        header = True
        for line in sample_file:
            if header:
                header = False
                continue
            sample = line.split()[0]
            samples.add(sample)
    return samples

unfiltered_samples = \
    load_samples(f'{ukb}/sample_qc/runs/{run_name}/combined.sample')
filtered_samples = \
    load_samples(f'{ukb}/sample_qc/runs/{run_name}/combined_unrelated.sample')

# from sample_id to sample_id
neighbors: Dict[str, Set[str]] = {}

for sample in unfiltered_samples:
    neighbors[sample] = set()

# Check if there is any relatedness
relatedness_fname = f'{ukb}/misc_data/ukbgene/ukb46122_rel_s488282.dat'
with open(relatedness_fname) as kinship_file:
    first = True
    for line in kinship_file:
        if first:
            first = False
            continue

        sample1, sample2 = line.split()[0:2]
        if sample1 in filtered_samples and sample2 in filtered_samples:
            error(f"Found related pair in the file {sample1}, {sample2}")
        if sample1 in unfiltered_samples and sample2 in unfiltered_samples:
            neighbors[sample1].add(sample2)
            neighbors[sample2].add(sample1)

# Check for maximality
for sample in unfiltered_samples:
    if sample in filtered_samples:
        continue
    found_error = True
    for neighbor in neighbors[sample]:
        if neighbor in filtered_samples:
            found_error = False
            break
    if found_error:
        error(f"Found a sample ({sample}) that was filtered but is not "
              f"related to anyone left")

# check for subset
for sample in filtered_samples:
    if sample not in unfiltered_samples:
        error(f"Found a sample ({sample}) in the unrelated set that's not "
              f"in the original set")

print("Success! No relatedness found. Remaining unrelated sample set is a "
      "maximal subset of the original.")
