#!/usr/bin/env python3

import collections
import os
import sqlite3
import sys
import time

import cyvcf2

import python_utils as utils

ukb = os.environ['UKB']

print('Loading refs ... ', end ='', flush=True)
start = time.time()
refs = utils.load_reference(f'{ukb}/../../resources/dbase/human/hg19/hg19.fa')

def get_next_strs(next_var, vcf):
    if next_var is None:
        for var in vcf:
            if var.ID not in str_ids:
                continue
            next_var = var
            break

    pos = next_var.INFO['START']
    border = max(min_border, next_var.INFO['PERIOD']*period_mult)
    curr_STRs = [(
        pos,
        next_var.REF[(pos-next_var.POS):(next_var.INFO['END']-next_var.POS+1)],
        next_var.INFO['PERIOD']
    )]
    curr_range = range(pos - border, next_var.INFO['END'] + 1 + border)
    del next_var

    for var in vcf:
        if var.ID not in str_ids:
            continue
        # filter duplicate STRs
        if var.INFO['START'] == curr_STRs[-1][0]:
            continue
        pos = var.INFO['START']
        border = max(min_border, var.INFO['PERIOD']*period_mult)
        if pos - border >= curr_range.stop:
            return (curr_STRs, curr_range, var)
        else:
            curr_STRs.append((
                pos,
                var.REF[(pos-var.POS):(var.INFO['END']-var.POS+1)],
                var.INFO['PERIOD']
            ))
            curr_range = range(curr_range.start, var.INFO['END'] + 1 + border)

    return (curr_STRs, curr_range, None)

def standardize(kmer):
    options = set()
    for i in range(len(kmer)):
        options.add(kmer[i:] + kmer[:i])
    return min(options)

def infer_repeat_unit(seq, period, pos):
    kmer_counts = collections.Counter()
    for i in range(len(seq) - period + 1):
        kmer = standardize(seq[i:(i+period)])
        kmer_counts[kmer] += 1
    best_kmers = kmer_counts.most_common(2)
    if len(best_kmers) == 1:
        return best_kmers[0][0]
    (top_kmer, top_count), (next_kmer, next_count) = best_kmers
    if next_count*2 >= top_count:
        return None
    return top_kmer

def is_pure_repeats(seq, repeat_unit):
    slen = len(seq)
    rlen = len(repeat_unit)
    for start_idx in range(rlen):
        pure = True
        for seq_idx in range(slen):
            if seq[seq_idx] != repeat_unit[(seq_idx + start_idx) % rlen]:
                pure = False
                break
        if pure:
            return True
    return False

for chrom in range(1, 23):
    ref_chrom = refs[chrom].upper()
    print(f'done. Time {time.time() - start}sec', flush=True)

    with open(f'{ukb}/snpstr/str_ids.txt') as str_ids_file:
        lines = str_ids_file.readlines()
    str_ids = set(line.strip() for line in lines)

    vcf = cyvcf2.VCF(f'{ukb}/snpstr/info_field/chr{chrom}.vcf.gz')
    snp_itr = sqlite3.connect(
        f"{ukb}/array_imputed/ukb_imp_chr{chrom}_v3.bgen.bgi"
    ).execute(
        "SELECT position, allele1, allele2 FROM Variant"
    )

    min_border = 3
    period_mult = 2


    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # don't handle reverse complement because all SNPs
    # are written according to the forward strand

    curr_STRs, curr_range, next_cyvcf2_str_var = get_next_strs(None, vcf)

    with open(f'{ukb}/finemapping/str_imp_snp_overlaps/chr{chrom}_to_examine.tab', 'w') as to_examine, \
            open(f'{ukb}/finemapping/str_imp_snp_overlaps/chr{chrom}_to_filter.tab', 'w') as to_filter, \
            open(f'{ukb}/finemapping/str_imp_snp_overlaps/chr{chrom}_to_keep.tab', 'w') as to_keep:
        header = 'ref_start\tref\t[str_pos, repeat_unit/period]\tsnp_pos\tsnp_ref\tsnp_alt\n'
        to_examine.write(header)
        to_filter.write(header)
        to_keep.write(header)
        chrom_done = False
        for pos, a1, a2 in snp_itr:
            if len(a1) == 0 or len(a2) == 0:
                raise ValueError("Don't want zero len alleles")
            while curr_range.stop <= pos:
                if next_cyvcf2_str_var is None:
                    chrom_done = True
                    break
                curr_STRs, curr_range, next_cyvcf2_str_var = get_next_strs(next_cyvcf2_str_var, vcf)
            if chrom_done:
                break
            end_pos_incl = pos + len(a1) - 1
            if curr_range.start > end_pos_incl:
                continue

            # there is an overlap, so decide: keep, filter or examine

            # values I will print
            print_STRs = []
            for STR in curr_STRs:
                repeat_unit = infer_repeat_unit(STR[1], STR[2], STR[0])
                if repeat_unit is not None:
                    print_STRs.append((STR[0], repeat_unit))
                else:
                    print_STRs.append((STR[0], f'period:{STR[2]}'))
            ref_seq = ''
            for idx in range(min(curr_range.start, pos), max(curr_range.stop, pos + len(a1))):
                if pos <= idx <= end_pos_incl:
                    ref_seq += ref_chrom[idx - 1].lower()
                else:
                    ref_seq += ref_chrom[idx - 1]
            out_str = f'{curr_range.start}\t{ref_seq}\t{print_STRs}\t{pos}\t{a1}\t{a2}\n'

            if len(curr_STRs) > 1:
                to_examine.write(out_str)
                continue

            curr_STR = curr_STRs[0]
            # no change in # bp can't be simple STR expansions or contractions
            # (unless we're swapping from one STR to another, but there's only
            # one STR in this region) so don't filter it
            if len(a1) == len(a2):
                to_keep.write(out_str)
                continue

            # since this is a single STR, if the SNP is not in the bounds of the STR, keep it
            if (pos + len(a1) - 1 < curr_STR[0] - 1 or
                pos > curr_STR[0] + len(curr_STR[1])):
                to_keep.write(out_str)
                continue

            if repeat_unit is None:
                to_examine.write(out_str)
                continue

            # keep any alternate alleles which obviously insert or delete impurities
            obvious_impurity = False
            for char in a2:
                if char not in repeat_unit and char not in a1:
                    obvious_impurity = True
                    break
            for char in a1:
                if char not in repeat_unit and char not in a2:
                    obvious_impurity = True
                    break
            if obvious_impurity:
                to_keep.write(out_str)
                continue

            # deletions of pure repeats in the STR are contractions and should be filtered
            # check for deletions where the first or the last characters match (for ease)
            # they should be in the STR proper
            if len(a2) == 1:
                deleted = None
                if a1[0] == a2 and pos >= curr_STR[0] -1:
                    deleted = a1[1:]
                elif a1[-1] == a2 and pos + len(a1) -1 <= curr_STR[0] + len(curr_STR[1]):
                    deleted = a1[:-1]
                if deleted is not None and is_pure_repeats(deleted, repeat_unit):
                    # check if the deleted portion is pure repeat (don't worry
                    # about correct rotation tho)
                    to_filter.write(out_str)
                    continue

            # insertions of pure repeats in the STR are expansions and should be filtered
            # check for insertions where the first or the last characters match (for ease)
            # they should be in the STR proper
            if len(a1) == 1:
                if a2[0] == a1 and pos >= curr_STR[0] -1:
                    inserted = a2[1:]
                elif a2[-1] == a1 and pos <= curr_STR[0] + len(curr_STR[1]):
                    inserted = a2[:-1]
                if inserted is not None and is_pure_repeats(inserted, repeat_unit):
                    # check if the inserted portion is pure repeat (don't worry
                    # about correct rotation tho)
                    to_filter.write(out_str)
                    continue

            to_examine.write(out_str)

            # TODO: maybe don't filter partial contraction or insertion in the middle?
            # either full contraction/insertion in the middle or partial at the edges?
