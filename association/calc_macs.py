#!/usr/bin/env python3

import argparse
import subprocess as sp
import sys
import tempfile

def imputed_snp_mac(outdir, plink_command, pfile_prefix, chrom, start, end, sample_fname):
    out = sp.run(
        f"cd {outdir} && "
        f"{plink_command} "
        f"--pfile {pfile_prefix} "
        f"--chr {chrom} "
        f"--from-bp {start} "
        f"--to-bp {end} "
        f"--keep {sample_fname} "
        '--memory 7900 '
        '--threads 4 '
        "--freq counts cols=+pos",
        shell=True,
        capture_output=True
    )

    print('stdout')
    print(out.stdout.decode(), flush=True)
    print('stderr', file=sys.stderr)
    print(out.stderr.decode(), file=sys.stderr, flush=True)
    out.check_returncode()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('plink_command')
    parser.add_argument('pfile_prefix')
    parser.add_argument('chrom')
    parser.add_argument('start', type=int)
    parser.add_argument('end', type=int)
    parser.add_argument('sample_file')
    args = parser.parse_args()

    with tempfile.NamedTemporaryFile('w') as temp_samps:
        with open(args.sample_file) as samps:
            # manually handle header
            temp_samps.write('FID ID')
            next(samps)
            for line in samps:
                line = line.strip()
                temp_samps.write(f'{line} {line}\n')
        
        imputed_snp_mac(args.outdir, args.plink_command, args.pfile_prefix, args.chrom, args.start, args.end, temp_samps.name)

if __name__ == '__main__':
    main()
