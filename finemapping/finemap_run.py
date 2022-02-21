#!/usr/bin/env python3

import argparse
import os
import subprocess as sp
import sys
import time

import h5py

import python_file_utils as file_utils

ukb = os.environ['UKB']

def run_finemap(outdir, *, prob_conv_sss_tol, prior_std, prior_snps):
    out = sp.run(
        f'{ukb}/utilities/finemap/finemap_v1.4_x86_64 --sss '
        f'--in-files {outdir}/finemap_input.master '
        '--log '
        '--n-configs-top 100 '
        '--n-threads 2 '
        #'--n-causal-snps 20 '
        '--n-causal-snps 40 ' # TODO undo this
        f'--prob-conv-sss-tol {prob_conv_sss_tol} '
        f'--prior-std {prior_std} ' +
        ('' if not prior_snps else '--prior-snps '),
        shell=True,
        capture_output=True
    )

    print(out.stdout.decode(), flush=True)
    print(out.stderr.decode(), file=sys.stderr, flush=True)
    out.check_returncode()

    time.sleep(10)

    if os.stat(f'{outdir}/finemap_output.snp').st_size == 0:
        raise RuntimeError("FINEMAP did not produce output!")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('--prob-conv-sss-tol', default=0.001, type=float)
    parser.add_argument('--prior-std', default=0.05, type=float)
    parser.add_argument('--prior-snps', default=False, action='store_true')
    args = parser.parse_args()

    run_finemap(
        args.outdir,
        prob_conv_sss_tol=args.prob_conv_sss_tol,
        prior_std=args.prior_std,
        prior_snps=args.prior_snps
    )

if __name__ == '__main__':
    main()
