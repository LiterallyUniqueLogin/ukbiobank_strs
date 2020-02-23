import argparse
import functools
import glob
import os

# files must already exist to be listed

UKB = os.environ['UKB']


def list_of_files(run_name, chrom, pos=False):
    file_type = 'samples'
    if pos:
        file_type = 'pos'
    files = glob.glob(f"{UKB}/str_imputed/runs/{run_name}/"
                      f"batches/chr{chrom}_{file_type}_*.vcf.gz")

    def num_string_comparator(a, b):
        if len(a) != len(b):
            return len(a) - len(b)
        elif a < b:
            return -1
        elif a == b:
            return 0
        else:
            return 1
    return sorted(files, key=functools.cmp_to_key(num_string_comparator))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_name",
                        help="will look for the vcfs to merge in "
                             "{UKB}/str_imputed/runs/{run_name}/"
                             "vcf_batches/chr{chr}_samples_*.vcf.gz")
    parser.add_argument("chr", help="the chromosome number")
    parser.add_argument("--pos", action="store_true")

    args = parser.parse_args()

    files = list_of_files(args.run_name, args.chr, pos=args.pos)

    for file in files:
        print(file)


if __name__ == "__main__":
    main()
