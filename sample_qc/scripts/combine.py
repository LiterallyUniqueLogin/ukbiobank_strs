#!/usr/bin/env python3

"""Combine sample file specifications into a single sample file."""
import argparse
import datetime
import sys

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('outsamples')
parser.add_argument('outreadme')
parser.add_argument('start_samples_fname')
parser.add_argument('remove_fnames', nargs='+')
parser.add_argument('--pheno-file')
args = parser.parse_args()
pheno_file = args.pheno_file

def error(msg):
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(-1)

with open(args.outreadme, 'w') as readme:
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    readme.write(f"Run date: {today}\n")

    if pheno_file:
        readme.write(f"Loading list of samples with phenotyeps from {pheno_file}\n")
    else:
        readme.write('Not subsetting to samples for any specific phenotype')

    remove_fnames = args.remove_fnames
    assert len(remove_fnames) > 0

    readme.write(f"Starting sample file: {args.start_samples_fname}\n")
    readme.write("Remove fnames:\n")
    for remove_fname in remove_fnames:
        readme.write(remove_fname + '\n')

if pheno_file:
    phen_array = np.load(pheno_file)
    samples = set(phen_array[:, 0].astype(int).astype(str))

start_lines = {"ID_1 ID_2 missing sex", "ID", "ID_1"}
with open(args.start_samples_fname) as start_samples:
    header = True
    current_keep = {}
    for line in start_samples:
        line = line.strip()
        if line == "":
            continue
        if header:
            if line not in start_lines:
                error(
                    f"Expected the first line of keep file {args.start_samples_fname} to be "
                    f"one of the header lines {start_lines}, instead see {line}"
                )
            header = False
            continue
        sample_id = line.split()[0]
        current_keep[sample_id] = line

    if not pheno_file:
        samples = set(current_keep)
    else:
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

with open(args.outsamples, 'w') as out_file:
    out_file.write('ID\n')
    for sample in samples:
        out_file.write(sample + "\n")
