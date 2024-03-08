#!/bin/env python3

"""
Read VCFs of STRs or SNPs,
apply filters to them for quality control
and then return yield them
"""

import os
from typing import Optional, Set, Tuple
import sys

import bgen_reader
import cyvcf2
import numpy as np
import scipy

import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

allele_len_precision = 2
dosage_precision = 2
r2_precision = 2

def dict_str(d):
    out = '{'
    first = True
    for key in sorted(d.keys()):
        if not first:
            out +=', '
        first = False
        # make sure keys are quoted so that the resulting string
        # is valid JSON and can be parsed as a dictionary in javascript
        # using JSON.parse
        out += f'{repr(str(key))}: {repr(d[key])}'
    out += '}'
    return out.replace("'", '"').replace('(', '[').replace(')', ']').replace('nan', '"NaN"')

def clean_len_alleles(d):
    new_d = {}
    for key, val in d.items():
        new_key = round(key, allele_len_precision)
        if new_key not in new_d:
            new_d[new_key] = val
        else:
            new_d[new_key] += val
    return new_d

def clean_len_allele_pairs(d):
    new_d = {}
    for (k1, k2), val in d.items():
        new_key = (round(k1, allele_len_precision), round(k2, allele_len_precision))
        if new_key not in new_d:
            new_d[new_key] = val
        else:
            new_d[new_key] += val
    return new_d

def round_vals(d, precision):
    return {key: round(val, precision) for key, val in d.items()}

def load_strs(
    vcf_fname: str,
    region: str,
    samples: np.ndarray,
    details: bool = True,
    var_subset: Optional[Set[int]] = None,
    best_guesses = False,
    both_poses = False,
    ignore_overlap_start = True
):
    """
    Iterate over a region returning genotypes at STR loci.

    First yield is a tuple of names of the fields in details.
    Every subsequent yield is described in the yields section below.

    Parameters
    ----------
    vcf_fname:
    region:
        chr:start-end
    samples:
        A boolean array of length nsamples determining which samples are included
        (True) and which are not

    Yields
    ------
    dosages: Dict[float, np.ndarray]
        A dictionary from unique length alleles to 2D arrays of size (n_samples, 2)
        which contain the dosages of those alleles for each haplotype
        Length dosage are measured in number of repeats.

        None if locus_filtered is not None
    
        If best_guesses, then instead of that just an array nx2 of length alleles
    unique_alleles: np.ndarray
        Array of unique length alleles (measured in number of repeats),
        same length as the dosages dict
    chrom: str
        e.g. '13'
    pos: int
    locus_filtered:
        None if the locus is not filtered, otherwise
        a string explaining why.
        'MAC<20' if the minor allele dosage is less than 20
        after sample subsetting, per plink's standard.

        None if best_guesses.
    locus_details:
        tuple of strings with the same length as the first yield
        with the corresponding order.

        None if best_guesses.
        

    Notes
    -----
    Hardcalls mentioned in the locus details are phased best_guesses, and in some
    corner cases will not correspond to the maximum likelihood unphased allele.
    """

    if best_guesses:
        assert var_subset is not None and not details

    chrom, region_poses = region.split(':')
    region_start, _ = [int(pos) for pos in region_poses.split('-')]
    vcf = cyvcf2.VCF(vcf_fname)

    if details:
        yield (
            'motif',
            'period',
            'ref_len',
            'total_per_allele_dosages',
            'total_best_guess_alleles',
            'total_best_guess_genotypes',
            'subset_total_per_allele_dosages',
            'subset_total_best_guess_alleles',
            'subset_total_best_guess_genotypes',
            'subset_het',
            'subset_entropy',
            'subset_HWEP',
            'subset_allele_dosage_r2'
        )
    for record in vcf(region):
        if ignore_overlap_start and record.POS < region_start:
            # records that overlap this region but started before this region
            # should be considered part of the pervious region and not returned here
            continue
        if record.INFO.get('PERIOD') is None:
            # there are a few duplicate loci which I didn't identify for removal before doing the imputation
            # this identifies and removes the first duplication which has the same exact coordinate but is present twice
            continue
        if (record.CHROM, record.POS) in [('17', 80520458), ('1', 247747217), ('1', 247848392), ('21', 47741815), ('8', 145231731)]:
            # there are a few duplicate loci which I didn't identify for removal before doing the imputation
            # this identifies and removes the second and third duplications which have some difference in flanking bps but are otherwise the same
            continue

        # some records will be indexed with VCF pos, some will be with actual repeat start pos, either is fine
        if var_subset is not None and record.POS not in var_subset and record.INFO.get('START') not in var_subset:
            continue

        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype='hipstr')

        if both_poses:
            ret_pos = (trrecord.pos, record.POS)
        else:
            ret_pos = trrecord.pos

        len_alleles = [trrecord.ref_allele_length] + trrecord.alt_allele_lengths
        len_alleles = [round(allele_len, allele_len_precision) for allele_len in len_alleles]

        if best_guesses:
            yield (trrecord.GetLengthGenotypes()[samples, :-1], np.unique(len_alleles), trrecord.chrom, ret_pos, None, None)
            continue

        if details:
            total_dosages = {_len: 0 for _len in np.unique(len_alleles)}
            for p in (1, 2):
                ap = trrecord.format[f'AP{p}']
                total_dosages[len_alleles[0]] += np.sum(np.maximum(0, 1 - np.sum(ap, axis=1)))
                for i in range(ap.shape[1]):
                    total_dosages[len_alleles[i+1]] += np.sum(ap[:, i])

            # TODO this isn't right - these best guess calls are only best guess if there's only one
            # allele per length, doesn't take into account prob splitting over imperfections
            total_best_guess_alleles = clean_len_alleles(trrecord.GetAlleleCounts())
            total_best_guess_genotypes = clean_len_allele_pairs(trrecord.GetGenotypeCounts())

        if isinstance(samples, slice):
            assert samples == slice(None)
            n_subset_samples = trrecord.GetNumSamples()
        else:
            n_subset_samples = int(np.sum(samples))

        subset_allele_probs = {
            _len: np.zeros((n_subset_samples, 2)) for _len in np.unique(len_alleles)
        }

        for p in (1, 2):
            ap = trrecord.format[f'AP{p}']
            subset_allele_probs[len_alleles[0]][:, (p-1)] += \
                    np.maximum(0, 1 - np.sum(ap[samples, :], axis=1))
            for i in range(ap.shape[1]):
                subset_allele_probs[len_alleles[i+1]][:, (p-1)] += ap[samples, i]

        subset_total_dosages = {
            _len: np.sum(subset_allele_probs[_len]) for _len in subset_allele_probs
        }

        if details:
            subset_total_best_guess_alleles = clean_len_alleles(trrecord.GetAlleleCounts(samples))
            subset_total_best_guess_genotypes = clean_len_allele_pairs(trrecord.GetGenotypeCounts(samples))
            subset_best_guess_allele_freqs = clean_len_alleles(trrecord.GetAlleleFreqs(samples))

            subset_het = utils.GetHeterozygosity(subset_best_guess_allele_freqs)
            subset_entropy = utils.GetEntropy(subset_best_guess_allele_freqs)
            subset_hwep = utils.GetHardyWeinbergBinomialTest(
                subset_best_guess_allele_freqs,
                subset_total_best_guess_genotypes
            )

            # https://www.cell.com/ajhg/fulltext/S0002-9297(09)00012-3#app1
            # Browning, Brian L., and Sharon R. Browning. "A unified approach to genotype imputation and haplotype-phase inference for large data sets of trios and unrelated individuals." The American Journal of Human Genetics 84.2 (2009): 210-223.
            # appendix 1
            subset_allele_dosage_r2 = {}

            subset_best_guesses = np.around(trrecord.GetLengthGenotypes()[samples, :-1], allele_len_precision)
            for length in len_alleles:
                # calculate allele dosage r**2 for this length
                if length in subset_allele_dosage_r2:
                    continue

                calls = subset_best_guesses == length

                subset_allele_dosage_r2[length] = np.corrcoef(
                    calls.reshape(-1), subset_allele_probs[length].reshape(-1)
                )[0,1]**2

            locus_details = (
                trrecord.motif,
                str(len(trrecord.motif)),
                str(round(trrecord.ref_allele_length, allele_len_precision)),
                dict_str(round_vals(total_dosages, dosage_precision)),
                dict_str(total_best_guess_alleles),
                dict_str(total_best_guess_genotypes),
                dict_str(round_vals(subset_total_dosages, dosage_precision)),
                dict_str(subset_total_best_guess_alleles),
                dict_str(subset_total_best_guess_genotypes),
                str(subset_het),
                str(subset_entropy),
                str(subset_hwep),
                dict_str(round_vals(subset_allele_dosage_r2, r2_precision))
            )
        else:
            locus_details = None

        mac = list(subset_total_dosages.values())
        mac.pop(np.argmax(mac))

        if np.sum(mac) < 20:
            yield (
                None,
                np.unique(len_alleles),
                trrecord.chrom,
                ret_pos,
                'MAC<20',
                locus_details
            )
            continue

        yield (
            subset_allele_probs,
            np.unique(len_alleles),
            trrecord.chrom,
            ret_pos,
            None,
            locus_details
        )

def load_imputed_snps(bgen_fname: str,
                      mfi_fname: Optional[str],
                      region: str,
                      samples: np.ndarray,
                      info_thresh: Optional[float] = None,
                      apply_filter: bool = True,
                      details: bool = True,
                      var_subset: Optional[Set[Tuple[int, str, str]]] = None,
                      best_guesses=False): # pos, ref, alt
    """
    Iterate over a region returning genotypes at SNP loci.

    First yield is a tuple of names of the fields in details.
    Every subsequent yield is described in the yields section below.

    Parameters
    ----------
    bgen_fname:
    mfi_fname:
        if details is None False and  info_thresh is None, can be omitted
        in that case, will occassionally include snps which are constant
    region:
        chr:start-end
    samples:
        A boolean array of length nsamples determining which samples are included
        (True) and which are not
    info_thresh:
        Loci must have INFO scores of at least this thresh or will be filtered
    apply_filter:
        should the mac 20 cutoff be applied?

    Yields
    ------
    dosages: np.ndarray
        A 2D array of size (n_samples, 3)
        which contain the probabilites of each genotype.

        None if locus_filtered is not None

        If best_guesses, then instead of that just an array nx2 of length alleles
    unique_alleles: np.ndarray
        Array of unique length alleles, same length as the dosages dict
    chrom: str
        e.g. '13'
    pos: int
    locus_filtered:
        None if the locus is not filtered, otherwise
        a string explaining why.
        'MAF=0' if there is only one allele present
        in the calls at this locus.
        'MAC<20' if the minor allele dosage is less than 20
        after sample subsetting, per plink's standard
        'info<{thresh}' if an info thresh
        is set and this locus is under that

        None if best_guesses
    locus_details:
        tuple of strings with the same length as the first yield
        with the corresponding order.

        None if best_guesses

    Notes
    -----
    1st allele is always the reference for imputed SNPs, alt allele
    is always nonreference.
    (This is different than microarray SNPs!)

    (Also note that these genotypes are not phased so the ordering
    of the genotypes for heterozygous samples is arbitrary. This
    currently doesn't matter as hard calls are)
    """
    assert (not apply_filter) or info_thresh is None

    if best_guesses:
        assert var_subset is not None and not details

    chrom, posses = region.split(':')
    start, end = tuple(int(val) for val in posses.split('-'))
    bgen = bgen_reader.open_bgen(bgen_fname,
                                 verbose=False,
                                 allow_complex=True)

    if details:
        yield (
            'info',
            "total_per_allele_dosages",
            'total_best_guesses',
            'subset_total_per_allele_dosages',
            'subset_total_best_guesses',
            'subset_HWEP'
        )

    start_num = np.searchsorted(bgen.positions, start)

    try:
        if mfi_fname:
            mfi = open(mfi_fname)

        for variant_num_offset, pos in enumerate(bgen.positions[start_num:]):
            variant_num = variant_num_offset + start_num
            if details or info_thresh:
                try:
                    mfi_line = next(mfi)
                except StopIteration:
                    return
        
            if pos < start:
                continue
            if pos > end:
                break

            alleles = np.array(bgen.allele_ids[variant_num].split(','))
            if var_subset is not None and (pos, alleles[0], alleles[1]) not in var_subset:
                continue

            if best_guesses:
                probs, missing, ploidy = bgen.read(
                    variant_num,
                    return_missings=True,
                    return_ploidies=True
                )
                probs = np.squeeze(probs)[samples, :]
                gts = np.argmax(probs, axis=1)
                yield (gts, alleles, chrom, pos, None, None)
                continue

            if details or info_thresh:
                info_str = mfi_line.split()[-1]

                if info_str == 'NA':
                    yield (
                        None,
                        alleles,
                        chrom,
                        pos,
                        'MAF=0',
                        ('NA', '', '', '', '', '')
                    )
                    continue

            probs, missing, ploidy = bgen.read(
                variant_num,
                return_missings=True,
                return_ploidies=True
            )
            probs = np.squeeze(probs)

            # make sure the record looks as expected
            '''
            assert len(probs.shape) == 2
            assert probs.shape[1] == 3
            assert not bgen.phased[variant_num]
            assert np.all(ploidy == 2)
            assert not np.any(missing)
            '''

            dosages = probs[:, 1] + 2*probs[:, 2]
            subset_dosages = dosages[samples]
            subset_total_alt_dosage = np.sum(subset_dosages)
            subset_total_ref_dosage = \
                    2*subset_dosages.shape[0] - subset_total_alt_dosage

            if details:
                total_alt_dosage = np.sum(dosages)
                total_ref_dosage = 2*dosages.shape[0] - total_alt_dosage

                best_guesses = probs.argmax(axis=1)
                n_hom_ref = np.sum(best_guesses == 0)
                n_het = np.sum(best_guesses == 1)
                n_hom_alt = np.sum(best_guesses == 2)

                subset_best_guesses = best_guesses[samples]
                subset_n_hom_ref = np.sum(subset_best_guesses == 0)
                subset_n_het = np.sum(subset_best_guesses == 1)
                subset_n_hom_alt = np.sum(subset_best_guesses == 2)

                subset_n_ref_alleles = 2*subset_n_hom_ref + subset_n_het
                subset_frac_ref_alleles = (subset_n_ref_alleles)/(2*subset_best_guesses.shape[0])
                subset_exp_hom_frac = subset_frac_ref_alleles**2 + (1-subset_frac_ref_alleles)**2
                subset_hwep = scipy.stats.binom_test(subset_n_hom_ref + subset_n_hom_alt,
                                                     n=subset_best_guesses.shape[0],
                                                     p=subset_exp_hom_frac)

                locus_details = (
                    info_str,
                    "{" f"'ref': {total_ref_dosage}, 'alt': {total_alt_dosage}" "}",
                    '{' f"'hom_ref': {n_hom_ref}, 'het': {n_het}, 'hom_alt': {n_hom_alt}" '}',
                    '{' f"'ref': {subset_total_ref_dosage}, 'alt': {subset_total_alt_dosage}" '}',
                    '{' f"'hom_ref': {subset_n_hom_ref}, 'het': {subset_n_het}, 'hom_alt': {subset_n_hom_alt}" '}',
                    f'{subset_hwep}'
                )
            else:
                locus_details = None

            if apply_filter and (subset_total_alt_dosage < 20 or subset_total_ref_dosage < 20):
                yield (
                    None,
                    alleles,
                    chrom,
                    pos,
                    'MAC<20',
                    locus_details
                )
                continue

            if info_thresh is not None and float(info_str) < info_thresh:
                yield (
                    None,
                    alleles,
                    chrom,
                    pos,
                    f'info<{info_thresh}',
                    locus_details
                )
                continue

            yield (
                probs[samples, :],
                alleles,
                chrom,
                pos,
                None,
                locus_details
            )
    finally:
        if mfi_fname:
            mfi.close()

