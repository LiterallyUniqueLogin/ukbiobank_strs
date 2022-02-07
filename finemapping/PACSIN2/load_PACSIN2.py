#!/usr/bin/env python3

import glob
import os

import numpy as np

import load_and_filter_genotypes as lfg
from PACSIN2_varnames import varnames

ukb = os.environ['UKB']

reflens = [53,3,18,14,18]

def get_gt_itr(samples):
    yield ('motif', 'period', 'ref_len', 'subset_total_per_allele_dosages')
    for varname, reflen in zip(varnames[:5], reflens):
        d = {}
        fnames = glob.glob(f'{ukb}/finemapping/PACSIN2/gts/{varname}*')
        for fname in fnames:
            len_ = int(fname.split('.')[0].split('_')[-1])
            d[len_] = np.load(fname)[samples, :]
        pos = int(varname.split('_')[1])
        motif = varname.split('_')[2]
        if motif == 'compound':
            period = 2
        elif  '|' not in motif:
            period = len(motif)
        else:
            period = 2
        subset_total_per_allele_dosages = {}
        for len_, gts in d.items():
            subset_total_per_allele_dosages[len_] = np.sum(gts)
        yield (d, np.array(sorted(d.keys())), 22, pos, None, (motif, str(period), str(reflen), lfg.dict_str(subset_total_per_allele_dosages)))

