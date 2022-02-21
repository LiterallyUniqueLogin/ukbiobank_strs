#!/usr/bin/env python3

import argparse
import ast
import os
import shutil
import subprocess as sp
import sys
import tempfile

import numpy as np

import sample_utils
import load_and_filter_genotypes

ukb = os.environ['UKB']

def imputed_snp_mac(outdir, chrom, start, end, sample_fname):
    out = sp.run(
        f"cd {outdir} ; "
        f"{ukb}/utilities/plink2 "
        f"--pfile {ukb}/array_imputed/pfile_converted/chr{chrom} "
        f"--chr {chrom} "
        f"--from-bp {start} "
        f"--to-bp {end} "
        f"--keep {sample_fname} "
        '--memory 3900 '
        "--freq counts cols=pos,ref,alt,altfreq ",
        shell=True,
        capture_output=True
    )

    print('stdout')
    print(out.stdout.decode(), flush=True)
    print('stderr', file=sys.stderr)
    print(out.stderr.decode(), file=sys.stderr, flush=True)
    out.check_returncode()

def str_mac(outdir, imputation_run_name, chrom, start, end, sample_fname):
    arr = sample_utils.samples_array_with_indicator(sample_fname)
    idx = ~np.isnan(arr[:, 1])
    itr = load_and_filter_genotypes.load_strs(imputation_run_name, f'{chrom}:{start}-{end}', idx)
    detail_names = next(itr)
    per_allele_dosages_idx = detail_names.index('subset_total_per_allele_dosages')
    with tempfile.NamedTemporaryFile('w') as temp_out:
        temp_out.write('pos\tmac\n')
        for _, _, _, pos, is_filtered, details in itr:
            if is_filtered is not None:
                continue
            acs = list(ast.literal_eval(details[per_allele_dosages_idx]).values())
            acs.pop(np.argmax(acs))
            mac = np.sum(acs)
            temp_out.write(f'{pos}\t{mac}\n')
        temp_out.flush()
        shutil.copy(temp_out.name, f'{outdir}/str_macs.tab')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('imputation_run_name')
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
        
        imputed_snp_mac(args.outdir, args.chrom, args.start, args.end, temp_samps.name)
        str_mac(args.outdir, args.imputation_run_name, args.chrom, args.start, args.end, args.sample_file)

if __name__ == '__main__':
    main()
