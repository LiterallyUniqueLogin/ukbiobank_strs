import collections

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
        
