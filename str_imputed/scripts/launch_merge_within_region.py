import argparse
import os
import subprocess as sp
import sys


def do_launch(run_name, chrom, regions, command_line=True):
    def error(msg):
        if command_line:
            print(msg, file=sys.stderr)
            exit(-1)
        else:
            return msg

    if "TMPDIR" not in os.environ:
        error("Didn't set TMPDIR environment variable")

    tmpdir = os.environ["TMPDIR"]
    ukb = os.environ["UKB"]

    if regions:
        for region in regions:
            if (int(region) - 1) % 5000000 != 0:
                error("Each region should be an index that is some "
                      "multiple of 5 million plus 1")

    if not regions:
        regions = range(1, 250000000, 5000000)

    for start_pos in regions:
        sp.run(f'sed -e "s/%RUN_NAME%/{run_name}/g" '
               f'-e "s/%CHROM%/{chrom}/g" '
               f'-e "s/%POS%/{start_pos}/g" '
               f'{ukb}/str_imputed/scripts/merge_within_region.pbs > '
               f'{tmpdir}/merge_within_region_{run_name}_{chrom}_{start_pos}'
               f'.pbs',
               shell=True, check=True)
        sp.run(f"qsub {tmpdir}/"
               f"merge_within_region_{run_name}_{chrom}_{start_pos}.pbs",
               shell=True, check=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_name",
                        help="Output files will be put in "
                        "$UKB/str_imputed/runs/run_name/batches")
    parser.add_argument("chromosome_number",
                        help="the number of the chromosome to regional merge")
    parser.add_argument("--regions",
                        nargs="+",
                        help="a list of regions to rerun (no other regions "
                             "will be run), e.g. --nargs 5000001 210000001 ")

    args = parser.parse_args()
    run_name = args.run_name
    chrom = args.chromosome_number
    do_launch(run_name, chrom, args.regions)

if __name__ == "__main__":
    main()
