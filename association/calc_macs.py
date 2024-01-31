#!/usr/bin/env python3

import argparse

import bgen_reader
import numpy as np

import sample_utils

def imputed_snp_mac(outfile, bgen_fname, start, end, samples):
    bgen = bgen_reader.open_bgen(bgen_fname, verbose=False)
    poses = (bgen.positions >= start) & (bgen.positions <= end)
    n_poses = np.sum(poses)
    n_samples = np.sum(samples)
    print(f"N vars: {n_poses}")
    with open(outfile, 'w') as out:
        out.write('chrom\tpos\tref\talt\tmac\tmaf\n')
        for count, idx in enumerate(np.where(poses)[0]):
            print(f'{count}/{n_poses}', end='\r')
            probs = bgen.read(idx).squeeze()[samples, :]
            mac = np.sum(probs[:, 1] + 2*probs[:, 2])
            if mac > n_samples:
                mac = n_samples*2 - mac
            out.write('\t'.join([
                bgen.chromosomes[idx],
                str(bgen.positions[idx]),
                *bgen.allele_ids[idx].split(','),
                str(mac),
                str(mac/(n_samples*2))
            ]) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile')
    parser.add_argument('bgen') # assume this is the correct chrom
    parser.add_argument('start', type=int)
    parser.add_argument('end', type=int)
    parser.add_argument('sample_file')
    parser.add_argument('all_samples_file')
    args = parser.parse_args()

    samples = sample_utils.get_samples_idx(args.all_samples_file, args.sample_file)
    imputed_snp_mac(args.outfile, args.bgen, args.start, args.end, samples)

if __name__ == '__main__':
    main()
