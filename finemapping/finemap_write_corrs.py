#!/usr/bin/env python3

import argparse

import h5py

import python_file_utils as file_utils

def write_corrs(workdir, outdir):
    '''
    write all_variants.ld
    run finemap
    write all finemap_output files
    '''

    with open(f'{workdir}/all_variants.ld', 'w') as ld_file, \
            h5py.File(f'{outdir}/lds.h5') as ldsh5:
        lds = ldsh5['ld']
        n_variants = lds.shape[0]
        for i in range(n_variants):
            ld_file.write(f'{lds[i, 0]:.10f}')
            for j in range(1, n_variants):
                ld_file.write(f' {lds[i, j]:.10f}')
            ld_file.write('\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    args = parser.parse_args()

    with file_utils.temp_dir(args.outdir) as tempdir:
        write_corrs(tempdir, args.outdir)
        file_utils.move_files(tempdir, args.outdir)

if __name__ == '__main__':
    main()
