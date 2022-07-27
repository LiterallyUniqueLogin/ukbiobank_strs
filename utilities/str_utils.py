import collections

import numpy as np

# repeat unit related
def standardize(kmer):
    options = set()
    for i in range(len(kmer)):
        options.add(kmer[i:] + kmer[:i])
    return min(options)

def infer_repeat_unit(seq, period):
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

complement_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
def canonicalize(seq):
    seq = seq.upper()
    choices = set()
    for i in range(len(seq)):
        cycle = seq[i:] + seq[:i]
        choices.add(cycle)
        reverse = ''
        for j in range(len(seq)):
            reverse += complement_dict[cycle[len(seq)-j-1]]
        choices.add(reverse)
    return min(choices)

# gt related
# in bps
# var - cyvcf2 variant
# samp_idx - 1d array of booleans of length n_samples
# these are length dosages
def str_dosage_gts(var, samp_idx):
    lens = [len(var.REF)] + [len(alt) for alt in var.ALT]
    gts = np.zeros(np.sum(samp_idx))
    gts += lens[0]*(
        np.maximum(0, 1 - np.sum(var.format('AP1'), axis=1)) +
        np.maximum(0, 1 - np.sum(var.format('AP2'), axis=1))
    )[samp_idx]
    for i in range(1, len(lens)):
        gts += lens[i]*(var.format('AP1') + var.format('AP2'))[samp_idx, i-1]
    return gts

# var, samp_idx - same as above
# alleles - a list of allele position (0 -ref, 1 - first alt, etc.) that contain the imperfection
# returns dosages of the number of imperfections per person
def imperfection_dosage_gts(var, samp_idx, alleles):
    gts = np.zeros(np.sum(samp_idx))
    for allele in alleles:
        if allele > 0:
            gts += (var.format('AP1') + var.format('AP2'))[samp_idx, allele-1]
        else:
            assert allele == 0
            gts += (
                np.maximum(0, 1 - np.sum(var.format('AP1'), axis=1)) +
                np.maximum(0, 1 - np.sum(var.format('AP2'), axis=1))
            )[samp_idx]
    return gts

