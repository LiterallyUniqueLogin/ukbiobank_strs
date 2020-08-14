import argparse
import os
import subprocess as sp

ukb = os.environ['UKB']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("previous_results_location")
    parser.add_argument("new_location")

    args = parser.parse_args()

    newname = args.new_location
    prevname = args.previous_results_location

    if os.path.exists(newname):
        print(f"File already exists at {newname}")
        exit(-1)

    with open(newname, 'w') as new, open(prevname) as prev:
        iprev = iter(prev)
        new.write(next(iprev)) # header
        for result_line in iprev:
            splits = result_line.split()
            chrom, pos, height_p, bilirubin_p = \
                    (float(splits[idx]) for idx in [0, 1, 4, 6])
            chrom = int(chrom)
            pos = int(pos)
            if height_p != 0 or bilirubin_p != 0:
                new.write(result_line)
                continue
            print(f"Fixing locus {chrom}:{pos}")
            cmd = (f'source ~/.bashrc ; '
                   f'conda activate ukb ; '
                   f'python {ukb}/association/simple_forloop.py height first_pass '
                   f'{chrom}_{pos}_rerun {chrom}:{pos}')
            out = sp.run(cmd, shell=True, capture_output=True, text=True)
            if out.returncode != 0:
                print(out.stdout)
                print(out.stderr)
                print("Call to simple_forloop.py failed with retcode "
                      f"{out.returncode}")
                exit(1)
            with open(f'{ukb}/association/runs/{chrom}_{pos}_rerun/'
                      'results.txt') as new_result:
                inew_result = iter(new_result)
                next(inew_result) # throw out header
                new.write(next(inew_result))


if __name__ == '__main__':
    main()
