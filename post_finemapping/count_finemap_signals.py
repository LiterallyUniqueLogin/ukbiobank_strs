#!/usr/bin/env python3

import glob
import os

ukb = os.environ['UKB']

total = 0
for count, file in enumerate(glob.glob(f'{ukb}/wdl_cache/finemapping/*_FINEMAP_first_pass_*cred*')):
    if count % 100 == 0:
        print(count)
    with open(file) as f:
        n_vars = int(file.split('.')[-1][4:])
        header = next(f)
        assert 'causal' in header, file
        prob = float(header.split()[8])
        total += prob*n_vars
print(total) 
