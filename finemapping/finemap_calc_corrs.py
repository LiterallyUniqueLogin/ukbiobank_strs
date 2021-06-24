#!/usr/bin/env python3

import argparse
import math
import os
import shutil

import h5py
import numpy as np

ukb = os.environ['UKB']

def correlate_chunks(gts, lds, n_variants, chunk_len, chunk_idx1, chunk_idx2):
    slice1 = slice(chunk_idx1*chunk_len, min((chunk_idx1+1)*chunk_len, n_variants))
    len_slice1 = slice1.stop - slice1.start
    slice2 = slice(chunk_idx2*chunk_len, min((chunk_idx2+1)*chunk_len, n_variants))
    gt1s = gts[slice1, :]
    gt2s = gts[slice2, :]
    corrs = np.corrcoef(gt1s, gt2s)
    lds[slice1, slice2] = corrs[:len_slice1, len_slice1:]
    lds[slice2, slice1] = corrs[len_slice1:, :len_slice1]
    print(f"Done with correlating chunks {chunk_idx1}, {chunk_idx2}", flush=True)
    assert not np.any(np.isnan(corrs))


def calc_corrs(workdir):
    '''
    write lds.h5 - dataset 'lds'
    '''

    chunk_len = 2**6

    with h5py.File(f'{workdir}/gts.h5') as gtsh5, \
            h5py.File(f'{workdir}/temp_lds.h5', 'w') as ldh5:
        gts = gtsh5['gts']
        n_variants = gts.shape[0]
        n_chunks = math.ceil(n_variants/chunk_len)

        lds = ldh5.create_dataset('ld', (n_variants, n_variants), chunks=True)
        for i in range(n_chunks):
            for j in range(i, n_chunks):
                correlate_chunks(gts, lds, n_variants, chunk_len, i, j)
    
    shutil.move(f'{workdir}/temp_lds.h5', f'{workdir}/lds.h5')

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

    calc_corrs(workdir)

if __name__ == '__main__':
    main()
