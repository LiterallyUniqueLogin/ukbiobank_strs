#!/usr/bin/env python3

import argparse
import math
import shutil

import h5py
import numpy as np

def correlate_chunks(gts, lds, n_variants, chunk_len, chunk_idx1, chunk_idx2):
    slice1 = slice(chunk_idx1*chunk_len, min((chunk_idx1+1)*chunk_len, n_variants))
    len_slice1 = slice1.stop - slice1.start
    slice2 = slice(chunk_idx2*chunk_len, min((chunk_idx2+1)*chunk_len, n_variants))
    gt1s = gts[:, slice1]
    gt2s = gts[:, slice2]
    corrs = np.corrcoef(gt1s.T, gt2s.T)
    lds[slice1, slice2] = corrs[:len_slice1, len_slice1:]
    lds[slice2, slice1] = corrs[len_slice1:, :len_slice1]
    print(f"Done with correlating chunks {chunk_idx1}, {chunk_idx2}", flush=True)
    assert not np.any(np.isnan(corrs))

def calc_corrs(outdir, gts_fname):
    '''
    write lds.h5 - dataset 'lds'
    '''

    chunk_len = 2**6

    with h5py.File(gts_fname) as gtsh5, \
            h5py.File(f'{outdir}/lds.h5.tmp', 'w') as ldh5:
        gts = gtsh5['gts']
        n_variants = gts.shape[1]
        n_chunks = math.ceil(n_variants/chunk_len)

        lds = ldh5.create_dataset('ld', (n_variants, n_variants), chunks=True)
        for i in range(n_chunks):
            for j in range(i, n_chunks):
                correlate_chunks(gts, lds, n_variants, chunk_len, i, j)

    shutil.move(f'{outdir}/lds.h5.tmp', f'{outdir}/lds.h5')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('gts')
    args = parser.parse_args()

    calc_corrs(args.outdir, args.gts)

if __name__ == '__main__':
    main()
