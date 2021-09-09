#!/usr/bin/env python3

import glob
import os

import numpy as np
import sortedcontainers

ukb = os.environ['UKB']

def get_pheno(fname):
    return '_'.join(fname.split('/')[-1].split('_')[:-2])

def load_peaks(fname):
    with open(fname) as pheno_file:
        next(pheno_file) # skip header
        peaks = []
        for line in pheno_file:
            split = line.split('\t')
            if split[5] == 'True':
                assert split[4] == 'True'
            if split[5] != '':
                assert float(split[3]) == 0
            peaks.append(split[2:6])
        # variant_type, pval, tagged_by_other, matched_p
        return peaks

def summarize_stat(peak_dict, name, predicate, percent_name=None, percent_predicate=None):
    assert bool(percent_name) == bool(percent_predicate)

    summary = {pheno: sum(predicate(variant_type, pval, tagged, matched_min_p) for variant_type, pval, tagged, matched_min_p in peaks) for pheno, peaks in peak_dict.items()}
    print(f'# {name}: {np.mean(list(summary.values())):.2f} ({np.std(list(summary.values())):.2f})')
    min_stat_pheno = min(summary, key=summary.get)
    min_stat = summary[min_stat_pheno]
    n_min = np.unique(list(summary.values()), return_counts=True)[1][0]
    if n_min > 1:
        min_str = f'{n_min} phenos'
    else:
        min_str = min_stat_pheno
    max_stat_pheno = max(summary, key=summary.get)
    max_stat = summary[max_stat_pheno]
    n_max = np.unique(list(summary.values()), return_counts=True)[1][-1]
    if n_max > 1:
        max_str = f'{n_max} phenos'
    else:
        max_str = max_stat_pheno
    print(f'min # {name}: {min_str}, {min_stat}   max: {max_str}, {max_stat}')

    if percent_name:
        denominator_summary =  {pheno: sum(percent_predicate(variant_type, pval, tagged, matched_min_p) for variant_type, pval, tagged, matched_min_p in peaks) for pheno, peaks in peak_dict.items()}
        perc_summary = {pheno: summary[pheno]/denominator_summary[pheno] for pheno in peak_dict}
        print(f'% {name} out of {percent_name}: {np.mean(list(perc_summary.values())) * 100:.2f}% ({np.std(list(perc_summary.values())) * 100:.2f}%)')
        min_stat_pheno = min(perc_summary, key=perc_summary.get)
        min_stat = perc_summary[min_stat_pheno]
        n_min = np.unique(list(perc_summary.values()), return_counts=True)[1][0]
        if n_min > 1:
            min_str = f'{n_min} phenos'
        else:
            min_str = min_stat_pheno
        max_stat_pheno = max(perc_summary, key=perc_summary.get)
        max_stat = perc_summary[max_stat_pheno]
        n_max = np.unique(list(perc_summary.values()), return_counts=True)[1][-1]
        if n_max > 1:
            max_str = f'{n_max} phenos'
        else:
            max_str = max_stat_pheno
        print(f'min % {name}: {min_str}, {min_stat*100:.2f}%   max: {max_str}, {max_stat*100:.2f}%')
        print(f'total: {sum(summary.values())/sum(denominator_summary.values())*100:.2f}%')

    print('')

def print_stats(peak_dict, group_name, display_largest_peaks=False):
    print(f'----------{group_name}----------')
    summarize_stat(
        peak_dict,
        'all peaks',
        lambda _, __, ___, ____ : True
    )
    for var in ('SNP', 'STR'):
        summarize_stat(
            peak_dict,
            f'tagged by {var} variant',
            lambda variant_type, _, tagged, __: variant_type == var or tagged == 'True',
            'all peaks',
            lambda _, __, ___, ____ : True
        )
        summarize_stat(
            peak_dict,
            f'{var} peaks',
            lambda variant_type, _, __, matched_min_p: variant_type == var and matched_min_p != 'True',
            'all peaks',
            lambda _, __, ___, ____ : True
        )
        summarize_stat(
            peak_dict,
            f'unmatched {var} peaks',
            lambda variant_type, _, tagged, matched_min_p: variant_type == var and tagged == 'False',
            f'{var} peaks',
            lambda variant_type, _, __, matched_min_p: variant_type == var and matched_min_p != 'True',
        )

    if display_largest_peaks:
        summarize_stat(
            peak_dict,
            'shared min_p peaks',
            lambda _, __, ___, matched_min_p: matched_min_p == 'True'
        )
        summarize_stat(
            peak_dict,
            'STR min_p peaks',
            lambda variant_type, __, ___, matched_min_p: variant_type == 'STR' and matched_min_p == 'False'
        )
    
def main():
    all_peaks = {}
    for fname in glob.glob(f'{ukb}/signals/peaks/*_500000_5e-8.tab'):
        all_peaks[get_pheno(fname)] = load_peaks(fname)
    n_phenos = len(all_peaks)
    print(f'# phenos: {n_phenos}')
    print_stats(all_peaks, 'regular', display_largest_peaks=True)

    str_peak_ranks_skip_shared_max_phenos = {}
    str_peak_ranks_ignore_shared_max = {}
    for pheno, peak_set in all_peaks.items():
        ranks = sortedcontainers.SortedSet((float(p_val), 'False' if matched_min_p == 'True' else 'True' if matched_min_p == 'False' else '', variant_type) for variant_type, p_val, _, matched_min_p in peak_set)
        shared_max = False
        rank = 0
        for p_val, not_matched_min_p, variant_type in ranks:
            if p_val == 0 and not_matched_min_p == 'False':
                shared_max = True
            else:
                rank += 1
                if variant_type == 'STR':
                    break
        str_peak_ranks_ignore_shared_max[pheno] = rank
        if not shared_max:
            str_peak_ranks_skip_shared_max_phenos[pheno] = rank
    assert len(str_peak_ranks_skip_shared_max_phenos) == 7
    assert len(str_peak_ranks_ignore_shared_max) == n_phenos

    print(str_peak_ranks_skip_shared_max_phenos)

    for data, name in (
        (str_peak_ranks_skip_shared_max_phenos, 'rank of best STR peak (phenos with shared max peaks skipped)'),
        (str_peak_ranks_ignore_shared_max, 'rank of best STR peak (ignore shared max peaks)')
    ):
        print(f'avg {name}: {np.mean(list(data.values())):.2f} ({np.std(list(data.values())):.2f})')
        min_stat_pheno = min(data, key=data.get)
        min_stat = data[min_stat_pheno]
        n_min = np.unique(list(data.values()), return_counts=True)[1][0]
        if n_min > 1:
            min_str = f'{n_min} phenos'
        else:
            min_str = min_stat_pheno
        max_stat_pheno = max(data, key=data.get)
        max_stat = data[max_stat_pheno]
        n_max = np.unique(list(data.values()), return_counts=True)[1][-1]
        if n_max > 1:
            max_str = f'{n_max} phenos'
        else:
            max_str = max_stat_pheno
        print(f'best rank: {min_str}, {min_stat}   worst rank: {max_str}, {max_stat}')

    print('')

    more_sig_peaks = {}
    for fname in glob.glob(f'{ukb}/signals/peaks/*_500000_5e-9.tab'):
        more_sig_peaks[get_pheno(fname)] = load_peaks(fname)
    assert n_phenos == len(more_sig_peaks)
    print_stats(more_sig_peaks, 'more sig')

    most_sig_peaks = {}
    for fname in glob.glob(f'{ukb}/signals/peaks/*_500000_5e-10.tab'):
        most_sig_peaks[get_pheno(fname)] = load_peaks(fname)
    assert n_phenos == len(most_sig_peaks)
    print_stats(most_sig_peaks, 'most sig')

    wider_peaks = {}
    for fname in glob.glob(f'{ukb}/signals/peaks/*_2000000_5e-8.tab'):
        wider_peaks[get_pheno(fname)] = load_peaks(fname)
    assert n_phenos == len(wider_peaks)
    print_stats(wider_peaks, 'wider')

if __name__ == '__main__':
    main()
