#!/usr/bin/env python3

import argparse
import shutil

import h5py

def write_corrs(outdir, lds_h5_fname):
    '''
    write all_variants.ld
    run finemap
    write all finemap_output files
    '''

    with open(f'{outdir}/all_variants.ld.tmp', 'w') as ld_file, \
            h5py.File(lds_h5_fname) as ldsh5:
        lds = ldsh5['ld']
        n_variants = lds.shape[0]
        for i in range(n_variants):
            ld_file.write(f'{lds[i, 0]:.10f}')
            for j in range(1, n_variants):
                ld_file.write(f' {lds[i, j]:.10f}')
            ld_file.write('\n')
    
    shutil.move(f'{outdir}/all_variants.ld.tmp', f'{outdir}/all_variants.ld')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('lds_h5')
    args = parser.parse_args()

    write_corrs(args.outdir, args.lds_h5)

if __name__ == '__main__':
    main()
