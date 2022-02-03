#!/usr/bin/env python3

import os
import re

ukb = os.environ['UKB']

with open(f'{ukb}/finemapping/PACSIN2/PACSIN2_alleles.txt') as f:
    ref, alt_alleles = next(f).strip().split()
    alt_alleles = alt_alleles.split(',')
    alleles = [ref] + alt_alleles

print('allele_num\tlen\tlen_A_reps\tlen_TA_reps\tlen_CA_reps\tlen_T(G|A)_reps\tmiddle_AAA\tmatch_num\n')
for allele_num, allele in enumerate(alleles):
    print(allele)
    m = re.match('^(A*)((?:TA)*)((?:CA)*)((?:T(?:G|A))*)$', allele)
    middle_a = False
    match_num=0
    if not m:
        m = re.match('^(A*)((?:TA)*)((?:CA)+GA(?:CA)+)((?:T(?:G|A))*)$', allele)
        match_num=1
    if not m:
        m = re.match('^(A*)((?:TA)+CA(?:TA)+)((?:CA)+)((?:T(?:G|A))*)$', allele)
        match_num=2
    if not m:
        middle_a = True
        m = re.match('^(A*)((?:TA)+)(AA(?:CA)+)((?:T(?:G|A))*TC)$', allele)
        match_num=3
    if not m:
        middle_a = False
        m = re.match('^(A*)((?:TA)*)((?:CA)+AA(?:CA)+)((?:T(?:G|A))*)$', allele)
        match_num=4
    if not m:
        m = re.match('^(A*)((?:TA)*)((?:CA)+)((?:T(?:G|A))+GA(?:T(?:G|A))+)$', allele)
        match_num=5
    assert m
    print(
        f'{allele_num}\t'
        f'{len(allele)}\t'
        f'{len(m[1])}\t'
        f'{len(m[2])}\t'
        f'{len(m[3])}\t'
        f'{len(m[4])}\t'
        f'{middle_a}\t'
        f'{match_num}'
    )

