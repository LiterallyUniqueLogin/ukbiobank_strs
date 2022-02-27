#!/usr/bin/env python3

import os

import matplotlib.pylot as plt
import pandas as pd
import upsetplot


ukb = os.environ['UKB']

upset_thresh = 0.8

df = pd.read_csv(f'{ukb}/post_finemapping/intermediate_results/finemapping_concordance.tab')

susie_tests = ['ratio', 'hardcall']
finemap_tests = ['ratio', 'conv_tol', 'total_prob', 'prior_std', 'mac']

for test in susie_tests:
    df[f'susie_alpha_{test}_indicator'] = df[f'susie_alpha_{test}']

for test in finemap_tests:
    df[f'finemap_pip_{test}_indicator'] = df[f'finemap_pip_{test}']

indicator_cols = (
    [f'susie_alpha_{test}_indicator' for test in susie_tests] +
    [f'finemap_pip_{test}_indicator' for test in finemap_tests]
)

upsetplot.from_indicators(indicator_cols, data = df)

plt.savefig(f'{ukb}/post_finemapping/results/all_upset.png')
