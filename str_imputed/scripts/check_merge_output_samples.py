import argparse
import glob
import gzip
import os
import sys

def collect_samples(ukb):
    samples = ""
    with open(glob.glob(f"{ukb}/microarray/*.sample")[0],
              'rt') as samples_file:
        samples_iter = iter(samples_file)
        next(samples_iter)
        first = True
        for line in samples_iter:
            line_head = line.split()[0]
            if line_head == "0":
                continue
            if first:
                first = False
            else:
                samples += "\t"
            samples += f"{line_head}_{line_head}"
    return samples


def do_check(run_name, chrom, command_line=True):
    ukb = os.environ['UKB']
    samples = collect_samples(ukb)
    failed_positions = set()
    for pos in range(1, 250000000, 5000000):
        end_pos = pos + 4999999
        file_name = (f"{ukb}/str_imputed/runs/{run_name}/batches/"
                     f"chr{chrom}_pos_{pos}_to_{end_pos}.vcf.gz")
        if command_line:
            print(f"Working on file {file_name}", end='\r')
        try:
            with gzip.open(file_name, 'rt') as imputed:
                imputed_iter = iter(imputed)
                line = '##'
                while line[:2] == '##':
                    line = next(imputed_iter)
                imputed_samples = line[46:].strip()
                if samples != imputed_samples:
                    failed_positions.add(pos)
                    if command_line:
                        sys.stdout.write("\033[K")
                        print(f"Samples in imputed file {file_name} not the "
                              "same as original callset samples!",
                              file=sys.stderr)
                        print("Original samples: ", samples, file=sys.stderr)
                        print("Imputed samples:",
                              imputed_samples,
                              file=sys.stderr)
                        sys.stdout.write("\033[K")
        except FileNotFoundError as exp:
            failed_positions.add(pos)
            if command_line:
                sys.stdout.write("\033[K")
                print(f"Couldn't find imputed file")
                print(exp)

    if command_line:
        if len(failed_positions) > 0:
            sys.stdout.write("\033[K")
            print(f"\nChromosome {chrom} failed sample check.",
                  file=sys.stderr)
            exit(-1)
        else:
            sys.stdout.write("\033[K")
            print(f"Samples are correct for all regional merge files for "
                  f"chromosome {chrom}")

    return failed_positions

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_name", help="The name of this run of imputation")
    parser.add_argument("chromosome_number",
                        help="The chromosome number to check")
    args = parser.parse_args()

    run_name = args.run_name
    chrom = args.chromosome_number
    do_check(run_name, chrom, command_line=True)

if __name__ == "__main__":
    main()
