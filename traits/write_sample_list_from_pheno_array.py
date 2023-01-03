#!/urs/bin/env python3

import argparse

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('pheno_data')
parser.add_argument('outfile')
args = parser.parse_args()

pheno_data = np.load(args.pheno_data)
with open(args.outfile, 'w') as sample_file:
    sample_file.write('ID\n')
    for ID in pheno_data[:, 0]:
        sample_file.write(str(int(ID)) + '\n')
