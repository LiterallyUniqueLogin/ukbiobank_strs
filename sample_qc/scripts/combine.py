#!/usr/bin/env python3

"""Combine sample file specifications into a single sample file."""
import argparse
import datetime
import glob
import os
import os.path
import sys

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('phenotype')
args = parser.parse_args()

ukb = os.environ['UKB']

rundir = f'{ukb}/sample_qc/runs/{args.phenotype}'
if not os.path.exists(rundir):
    os.mkdir(rundir)

def error(msg):
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(-1)

with open(f'{rundir}/README.txt', 'w') as readme:
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    readme.write(f"Run date: {today}\n")

    phen_sample_file = f'{ukb}/traits/phenotypes/{args.phenotype}.npy'
    readme.write(f"Loading list of samples with phenotyeps from {phen_sample_file}\n")

    keep_fnames = glob.glob(f'{ukb}/sample_qc/common_filters/keep/*sample')
    assert len(keep_fnames) > 0
    remove_fnames = glob.glob(f'{ukb}/sample_qc/common_filters/remove/*sample')
    assert len(remove_fnames) > 0

    readme.write("Keep fnames:\n")
    for keep_fname in keep_fnames:
        readme.write(keep_fname + '\n')
    readme.write("Remove fnames:\n")
    for remove_fname in remove_fnames:
        readme.write(remove_fname + '\n')

phen_array = np.load(phen_sample_file)
samples = set(phen_array[:, 0].astype(int).astype(str))

start_lines = {"ID_1 ID_2 missing sex", "ID", "ID_1"}
for keep_fname in keep_fnames:
    with open(keep_fname) as keep_file:
        header = True
        current_keep = {}
        for line in keep_file:
            line = line.strip()
            if line == "":
                continue
            if header:
                if line not in start_lines:
                    error(
                        f"Expected the first line of keep file {keep_fname} to be "
                        f"one of the header lines {start_lines}, instead see {line}"
                    )
                header = False
                continue
            sample_id = line.split()[0]
            current_keep[sample_id] = line

        samples = samples.intersection(current_keep)

for remove_fname in remove_fnames:
    with open(remove_fname) as remove_file:
        header = True
        for line in remove_file:
            line = line.strip()
            if line == "":
                continue
            if header:
                if line not in start_lines:
                    error(
                        f"Expected the first line of remove file {remove_fname} to be "
                        f"one of the header lines {start_lines}, instead see {line}"
                    )
                header = False
                continue
            sample = line.split()[0]
            if sample in samples:
                samples.remove(sample)

run_fname = f"{ukb}/sample_qc/runs/{args.phenotype}/combined.sample"
with open(run_fname, 'w') as out_file:
    out_file.write('ID\n')
    for sample in samples:
        out_file.write(sample + "\n")
