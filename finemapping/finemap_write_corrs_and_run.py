#!/usr/bin/env python3

import argparse
import os
import subprocess as sp
import sys

import h5py

ukb = os.environ['UKB']

def write_corrs_and_run(workdir):
    '''
    write all_variants.ld
    run finemap
    write all finemap_output files
    '''

    with open(f'{workdir}/all_variants.ld', 'w') as ld_file, \
            h5py.File(f'{workdir}/lds.h5') as ldsh5:
        lds = ldsh5['ld']
        n_variants = lds.shape[0]
        for i in range(n_variants):
            ld_file.write(f'{lds[i, 0]:.10f}')
            for j in range(1, n_variants):
                ld_file.write(f' {lds[i, j]:.10f}')
            ld_file.write('\n')

    out = sp.run(
        f'{ukb}/utilities/finemap/finemap_v1.4_x86_64 --sss '
        f'--in-files {workdir}/finemap_input.master '
        '--log '
        '--n-configs-top 100 '
        '--n-threads 2 '
        '--n-causal-snps 20',
        shell=True,
        capture_output=True
    )

    print(out.stdout.decode(), flush=True)
    print(out.stderr.decode(), file=sys.stderr, flush=True)
    out.check_returncode()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('chrom', type=int)
    parser.add_argument('start_pos', type=int)
    parser.add_argument('end_pos', type=int)
    args = parser.parse_args()

    phenotype = args.phenotype
    chrom = args.chrom
    start_pos = args.start_pos
    end_pos = args.end_pos
    assert start_pos < end_pos

    workdir = f'{ukb}/finemapping/finemap_results/{phenotype}/{chrom}_{start_pos}_{end_pos}'

    write_corrs_and_run(workdir)

if __name__ == '__main__':
    main()
