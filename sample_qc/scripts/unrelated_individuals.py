#!/usr/bin/env python3
"""
Find the unrelated sample subset for a run.

Writes the file
    $UKB/sample_qc/runs/{run_name}/combined_unrelated.sample
"""

import argparse
import datetime
import logging
import os
import subprocess as sp
import sys

parser = argparse.ArgumentParser()
parser.add_argument('run_name',
                    help='name of the run to look for combined.sample in')
args = parser.parse_args()
ukb = os.environ['UKB']

output_dir = f'{ukb}/sample_qc/runs/{args.run_name}/primus_output/'
os.makedirs(output_dir, exist_ok=True)
logging.basicConfig(filename=f"{output_dir}/script_log.log",
                    level=logging.DEBUG)

all_samples = set()
combined_floc = f'{ukb}/sample_qc/runs/{args.run_name}/combined.sample'
with open(combined_floc) as sample_file:
    header = True
    for line in sample_file:
        if header:
            header = False
            continue
        sample_id = line.strip()
        all_samples.add(sample_id)
logging.info(f"Done reading original sample file. Total num samples: "
             f"{len(all_samples)}")

# all sample numbers that appear in the related samples file
related_samples = set()

kinship_subset_floc = output_dir + "/kinship_subset.dat"
original_kinship_floc = (f'{ukb}/misc_data/ukbgene/ukb46122_rel_s488282.dat')
with open(original_kinship_floc) as original_kinship_file, \
        open(kinship_subset_floc, 'w') as kinship_subset_file:
    first = True
    for line in original_kinship_file:
        if first:
            first = False
            kinship_subset_file.write(line)
            continue

        sample1, sample2 = line.split()[0:2]
        if not (sample1 in all_samples and sample2 in all_samples):
            continue
        related_samples.add(sample1)
        related_samples.add(sample2)
        kinship_subset_file.write(line)

logging.info(f"Done creating subset kinship file. Num samples with relations: "
             f"{len(related_samples)}")

# FYI, The PRIMUS command produces an ungodly amount of auxiliary files when
# running - enough so that rm'ing them all takes a while
commandString = (
    f"{os.environ['SOURCE']}/PRIMUS_v1.9.0/bin/run_PRIMUS.pl "
    f"-i FILE={kinship_subset_floc} FID1=1 IID1=1 FID2=2 IID2=2 PI_HAT=5 "
    f"--no_PR -t 0.04419417382 -o {output_dir}/PRIMUS"
)
# the -t threshold equal to 1/2^(9/2) as described in the UKB nature paper

logging.info("PRIMUS command: " + commandString)

output = sp.run(commandString, shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                check=True)
logging.info("Done running PRIMUS")
logging.info("PRIMUS output: " + output.stdout.decode())
logging.info("PRIMUS error (if any): " + output.stderr.decode())

output_floc = (f"{ukb}/sample_qc/runs/{args.run_name}/"
               f"combined_unrelated.sample")
with open(output_floc, 'w') as output_file:
    output_file.write("ID\n")
    # write out all the samples that aren't related to anyone
    for sample in all_samples:
        if sample not in related_samples:
            output_file.write(sample + '\n')

    # write out all the samples that were selected among the related ones
    unrelated_floc = \
        f'{output_dir}/PRIMUS/kinship_subset.dat_maximum_independent_set'
    with open(unrelated_floc) as unrelated_file:
        first = True
        for line in unrelated_file:
            if first:
                first = False
                continue
            output_file.write(line.strip() + "\n")
