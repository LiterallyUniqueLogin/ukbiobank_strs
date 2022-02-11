#!/usr/bin/env python3
"""
Find the unrelated sample subset for a run.
"""

import argparse
import os
import shutil
import subprocess as sp
import tempfile

import numpy as np

def primus_command(infile, outdir):
    # FYI, The PRIMUS command produces an ungodly amount of auxiliary files when
    # running - enough so that rm'ing them all takes a while
    return (
        f"{os.environ['SOURCE']}/PRIMUS_v1.9.0/bin/run_PRIMUS.pl "
        f"-i FILE={infile} FID1=1 IID1=1 FID2=2 IID2=2 PI_HAT=5 "
        f"--no_PR -t 0.04419417382 -o {outdir}"
    )
    # the -t threshold equal to 1/2^(9/2) as described in the UKB nature paper

parser = argparse.ArgumentParser()
parser.add_argument('outfile')
parser.add_argument('kinship_file')
parser.add_argument('pheno_file')
parser.add_argument('--binary-pheno')
args = parser.parse_args()

with tempfile.TemporaryDirectory() as output_dir:
    print('Temp dir name: ', output_dir)

    all_samples = {str(int(id_)) for id_ in np.load(args.pheno_file)[:, 0]}

    print(f"Done reading original sample file. Total num samples: {len(all_samples)}")

    kinship_subset_fname = f'{output_dir}/kinship_subset.dat'

    if args.binary_pheno:
        # for case control data, prioritize cases over controls
        # and prioritize maximizing number of cases over maximizing
        # total number of cases and controls
        data = np.load(args.binary_pheno)
        # assert that cases are a low enough percentage of total samples
        # so that this type of prioritization makes sense
        assert np.sum(data[:, 1]) <= .2*data.shape[0]
        cases = set(str(int(ID)) for ID in data[data[:, 1] == 1, 0])
        print(f'N_cases: {len(cases)}')

        kinship_only_cases_fname = f'{output_dir}/cases_kinship_subset.dat'

        # collect all related cases
        any_lines = False
        related_cases = set()
        with open(args.kinship_file) as original_kinship_file, \
                open(kinship_only_cases_fname, 'w') as kinship_only_cases_file:
            kinship_only_cases_file.write(next(original_kinship_file))

            for line in original_kinship_file:
                sample1, sample2 = line.split()[0:2]
                if (
                    sample1 in all_samples and sample1 in cases and
                    sample2 in all_samples and sample2 in cases
                ):
                    kinship_only_cases_file.write(line)
                    related_cases.add(sample1)
                    related_cases.add(sample2)

        if not related_cases:
            print("No related cases")
        else:
            print(f"Done writing out cases_only relatedness file {kinship_only_cases_fname}.")
            print(f"N_related_cases: {len(related_cases)}")
            # get a maximal set of unrelated cases
            old_cases_primus_out = f'{output_dir}/case_PRIMUS_OLD'
            if os.path.exists(old_cases_primus_out):
                shutil.rmtree(old_cases_primus_out)
            cases_command_string = primus_command(kinship_only_cases_fname, f'{output_dir}/case_PRIMUS')

            output = sp.run(cases_command_string, shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                            check=True)
            print("Done running case PRIMUS")
            print("Case PRIMUS output: " + output.stdout.decode())
            print("Case PRIMUS error (if any): " + output.stderr.decode())

            unrelated_cases = set()
            unrelated_case_fname = \
                f'{output_dir}/case_PRIMUS/cases_kinship_subset.dat_maximum_independent_set'
            with open(unrelated_case_fname) as unrelated_case_file:
                next(unrelated_case_file) # skip header
                for line in unrelated_case_file:
                    unrelated_cases.add(line.strip())
            print(f'N_unrelated_cases: {len(unrelated_cases)}')

            # remove all the cases that aren't in the unrelated cases set
            all_samples = all_samples.difference(related_cases.difference(unrelated_cases))
            print(f'N controls + unrelated cases: {len(all_samples)}')

        # remove all the controls that are related to remaining cases
        with open(args.kinship_file) as original_kinship_file:
            next(original_kinship_file) # skip header

            for line in original_kinship_file:
                sample1, sample2 = line.split()[0:2]
                if   sample1 in all_samples and sample1 in cases and sample2 in all_samples:
                    all_samples.remove(sample2)
                elif sample2 in all_samples and sample2 in cases and sample1 in all_samples:
                    all_samples.remove(sample1)

        print(f'N unrelated cases and all controls unrelated to them: {len(all_samples)}')
        print("Done handling case/control specially.")

    # all sample numbers that appear in the related samples file
    related_samples = set()

    with open(args.kinship_file) as original_kinship_file, \
            open(kinship_subset_fname, 'w') as kinship_subset_file:
        kinship_subset_file.write(next(original_kinship_file)) # copy header

        for line in original_kinship_file:
            sample1, sample2 = line.split()[0:2]
            if not (sample1 in all_samples and sample2 in all_samples):
                continue
            related_samples.add(sample1)
            related_samples.add(sample2)
            kinship_subset_file.write(line)

    print(f"Done creating subset kinship file. Num samples with relations: {len(related_samples)}")

    primus_out_old = f'{output_dir}/PRIMUS_OLD'
    if os.path.exists(primus_out_old):
        shutil.rmtree(primus_out_old)
    command_string = primus_command(kinship_subset_fname, f'{output_dir}/PRIMUS')

    print("PRIMUS command: " + command_string)

    output = sp.run(command_string, shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                    check=True)
    print("Done running PRIMUS")
    print("PRIMUS output: " + output.stdout.decode())
    print("PRIMUS error (if any): " + output.stderr.decode())

    with open(args.outfile, 'w') as output_file:
        output_file.write("ID\n")
        # write out all the samples that aren't related to anyone
        for sample in all_samples:
            if sample not in related_samples:
                output_file.write(sample + '\n')

        # write out all the samples that were selected among the related ones
        unrelated_fname = \
            f'{output_dir}/PRIMUS/kinship_subset.dat_maximum_independent_set'
        with open(unrelated_fname) as unrelated_file:
            first = True
            for line in unrelated_file:
                if first:
                    first = False
                    continue
                output_file.write(line.strip() + "\n")
