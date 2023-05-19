#!/usr/bin/env python3

import argparse
import os
import subprocess as sp
import sys
import time

def run_finemap(outdir, finemap_command, *, n_causal_snps, prob_conv_sss_tol, prior_std, prior_snps):
    with open(f'{outdir}/finemap_input.z') as input_vars:
        n_vars = len(input_vars.readlines()) - 1
    causal_snps_flag = min(n_causal_snps, n_vars)

    out = sp.run(
        finemap_command +
        ' --sss '
        f'--in-files {outdir}/finemap_input.master '
        '--log '
        #'--n-configs-top 100 ' This was an error
        '--n-threads 2 '
        f'--n-causal-snps {causal_snps_flag} '
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
    parser.add_argument('finemap_command')
    parser.add_argument('--n-causal-snps', default=20, type=int)
    parser.add_argument('--prob-conv-sss-tol', default=0.001, type=float)
    parser.add_argument('--prior-std', default=0.05, type=float)
    parser.add_argument('--prior-snps', default=False, action='store_true')
    args = parser.parse_args()

    run_finemap(
        args.outdir,
        args.finemap_command,
        n_causal_snps=args.n_causal_snps,
        prob_conv_sss_tol=args.prob_conv_sss_tol,
        prior_std=args.prior_std,
        prior_snps=args.prior_snps
    )

if __name__ == '__main__':
    main()
