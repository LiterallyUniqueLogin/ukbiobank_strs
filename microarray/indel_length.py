#!/usr/bin/env python3

import collections
import sqlite3

import numpy as np
from statsmodels.stats.weightstats import DescrStatsW

# from https://enkre.net/cgi-bin/code/bgen/wiki/The%20bgenix%20index%20file%20format
biggest_deletion = 0
del_sizes = collections.defaultdict(int)
biggest_insertion = 0
ins_sizes = collections.defaultdict(int)
snp_count = 0

var_count = 0
for chrom in range(1, 23):
    index = sqlite3.connect(f"ukb_hap_chr{chrom}_v2.bgen.bgi")
    cursor = index.execute("select allele1, allele2, position from Variant")
    for a1, a2, pos in cursor:
        var_count += 1
        if var_count % 1000000 == 0:
            print(f'count: {var_count}')
            print(f'snp_count: {snp_count}')
            ins_stats = DescrStatsW(np.array(list(ins_sizes.keys())), weights=np.array(list(ins_sizes.values())))
            print(f'ins_count: {sum(list(ins_sizes.values()))}, ins_sizes: {ins_stats.mean:.4f} (+/-{ins_stats.std_mean:.4f}), biggest_ins: {biggest_insertion}')
            del_stats = DescrStatsW(np.array(list(del_sizes.keys())), weights=np.array(list(del_sizes.values())))
            print(f'del_count: {sum(list(del_sizes.values()))}, del_sizes: {del_stats.mean:.4f} (+/-{del_stats.std_mean:.4f}), biggest_del: {biggest_deletion}')

        if not set(a1.lower()).issubset('gcta') or not set(a2.lower()).issubset('gcta'):
            print(chrom, pos, a1, a2)
            assert False

        if len(a1) == len(a2) == 1:
            snp_count += 1
            continue

        oria1 = a1
        oria2 = a2
        while len(a1) > 0 and len(a2) > 0:
            if a1[0] == a2[0]:
                a1 = a1[1:]
                a2 = a2[1:]
            else:
                break
        while len(a1) > 0 and len(a2) > 0:
            if a1[-1] == a2[-1]:
                a1 = a1[:-1]
                a2 = a2[:-1]
            else:
                break

        if len(a1) == len(a2) == 1:
            snp_count += 1
            continue

        if not len(a1) == 0 and not len(a2) == 0:
            print(chrom, pos, oria1, oria2)
            assert False

        '''
        if len(a1) > 1 and len(a2) > 1 and not a1.startswith(a2) and not a2.startswith(a1):
            print(chrom, pos, a1, a2)
            assert False
        '''
        if len(a1) > len(a2):
            del_size = len(a1) - len(a2)
            del_sizes[del_size] += 1
            biggest_deletion = max(biggest_deletion, del_size)
        else:
            ins_size = len(a2) - len(a1)
            ins_sizes[ins_size] += 1
            biggest_insertion = max(biggest_insertion, ins_size)

print('--------------')
print(f'count: {var_count}')
print(f'snp_count: {snp_count}')
ins_stats = DescrStatsW(np.array(list(ins_sizes.keys())), weights=np.array(list(ins_sizes.values())))
print(f'ins_count: {sum(list(ins_sizes.values()))}, ins_sizes: {ins_stats.mean:.4f} (+/-{ins_stats.std_mean:.4f}), biggest_ins: {biggest_insertion}')
del_stats = DescrStatsW(np.array(list(del_sizes.keys())), weights=np.array(list(del_sizes.values())))
print(f'del_count: {sum(list(del_sizes.values()))}, del_sizes: {del_stats.mean:.4f} (+/-{del_stats.std_mean:.4f}), biggest_del: {biggest_deletion}')



