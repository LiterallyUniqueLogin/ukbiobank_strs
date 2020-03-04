# pylint: disable=C0103
"""Combine sample file specifications into a single sample file."""
import argparse
import os
import sys
from typing import Dict, List

parser = argparse.ArgumentParser()
parser.add_argument('--run_name')
parser.add_argument('-k', '--keep', action='append', nargs="+",
                    required=True,
                    help=("A list of sample files. Samples will be included"
                          " in the output only if they are present in all of "
                          "the keep files. At least one keep file must be "
                          "listed. This flag may be used more than once."))
parser.add_argument('-r', '--remove', action='append', nargs="+",
                    help=("A list of sample files. All the samples in these"
                          " files will be removed from the output."
                          " This flag may be used more than once."))

args = parser.parse_args()
keep_flocs: List[str] = sum(args.keep, [])
remove_flocs: List[str] = sum(args.remove, [])


def error(msg):
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(-1)


if len(keep_flocs) == 0:
    error(("Must specify at least one keep file. If you only wish to "
           "use remove files, simply specify the .sample file containing all "
           "the samples as the one keep file."))
    sys.exit(-1)

ukb = os.environ['UKB']


# Write out a log documenting which files we are pulling from
log_string = "Keep files:\n"
for floc in keep_flocs:
    log_string += f"{floc}\n"
log_string += "Remove files:\n"
for floc in remove_flocs:
    log_string += f"{floc}\n"

log_floc = f"{ukb}/sample_qc/runs/{args.run_name}/samples.log"
if os.path.exists(log_floc):
    with open(log_floc) as log_file:
        log = log_file.read()
        if log != log_string:
            error("Samples log already exists with different content!")
else:
    with open(log_floc, 'x') as log_file:
        log_file.write(log_string)

# Collect the samples
# Dict from ID to the sample line
# for that ID (e.g. which includes the gender)
sample_dict: Dict[str, str] = {}

first_keep = True
expected_first_line = "ID_1 ID_2 missing sex\n"
other_expected_starts = {"ID", "ID_1"}
for floc in keep_flocs:
    with open(floc) as keep_file:
        header = True
        current_keep = {}
        for line in keep_file:
            if line.strip() == "":
                continue
            sample_id = line.split()[0]
            if header:
                if first_keep and line != expected_first_line:
                    error(("Expected first line of the first keep file "
                           "'{floc}' to be the header line "
                           "'{expected_first_line}' instead of {line}"))
                elif sample_id not in other_expected_starts:
                    error(("Expected first line to be a header line that "
                           " begins with one of {other_expected_starts} in "
                           f"file {floc}, instead see {sample_id}"))
                header = False
                continue
            current_keep[sample_id] = line
        if first_keep:
            sample_dict = current_keep
            first_keep = False
        else:
            for sample_id in sample_dict:
                if sample_id not in current_keep:
                    del sample_dict[sample_id]

for floc in remove_flocs:
    with open(floc) as remove_file:
        header = True
        for line in remove_file:
            if line.strip() == "":
                continue
            sample_id = line.split()[0]
            if header:
                if sample_id not in other_expected_starts:
                    error(("Expected first line to be a header line that "
                           " begins with one of {other_expected_starts} in "
                           f"file {floc}, instead see {sample_id}"))
                header = False
                continue
            if sample_id in sample_dict:
                del sample_dict[sample_id]


run_floc = f"{ukb}/sample_qc/runs/{args.run_name}/combined.samples"
with open(run_floc, 'w') as out_file:
    out_file.write('ID_1 ID_2 missing sex\n')
    for sample_id in sample_dict:
        out_file.write(sample_dict[sample_id])
