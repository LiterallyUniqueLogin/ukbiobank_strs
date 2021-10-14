#!/usr/bin/env python3

import argparse
import bisect
import os
import pathlib

ukb = os.environ['UKB']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    args = parser.parse_args()
    phenotype = args.phenotype

    regions = []
    with open(f'{ukb}/signals/regions/{phenotype}.tab') as regions_itr:
        next(regions_itr)
        for line in regions_itr:
            regions.append(tuple(int(el) for el in line.strip().split('\t')))
    peaks = []
    with open(f'{ukb}/signals/peaks/{phenotype}_250000_5e-8.tab') as peaks_itr:
        next(peaks_itr)
        for line in peaks_itr:
            split = line.split('\t')
            peaks.append((int(split[0]), int(split[1]), 0))

    print('Checking each region contains a peak ...')
    for region in regions:
        if region[0] == 6 and region[1] <= 33.5e6 and region[2] >= 25e6:
            # ignore regions intersecting MHC
            continue
        idx = bisect.bisect(peaks, region)
        peak = peaks[idx]
        if (peak[0] != region[0] or
            peak[1] <= region[1] or
            peak[1] >= region[2]):
            print(peak, region)
            assert False

    print('Checking each peak is in a region ...')
    for peak in peaks:
        if peak[0] == 6 and 25e6 <= peak[1] <= 33.5e6:
            # ignore MHC peaks
            continue
        idx = bisect.bisect(regions, peak)
        region = regions[idx - 1]
        if (peak[0] != region[0] or
            peak[1] <= region[1] or
            peak[1] >= region[2]):
            print(peak, region)
            assert False

    print('Done.')

    pathlib.Path(f'{ukb}/signals/comparison/{phenotype}.done').touch()

if __name__ == '__main__':
    main()
