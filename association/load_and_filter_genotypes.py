#!/bin/env python3

"""
Read VCFs of STRs or SNPs,
apply filters to them for quality control
and then return yield them
"""

import os
import os.path
from typing import Optional, Set, Tuple
import sys

import bgen_reader
import cyvcf2
import numpy as np
import scipy

ukb = os.environ['UKB']

sys.path.insert(0, f'{ukb}/../trtools/repo')

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

def load_strs(imputation_run_name: str,
              region: str,
              samples: np.ndarray,
              details: bool = True,
              var_subset: Optional[Set[int]] = None,
              hardcalls = False):
    """
    Iterate over a region returning genotypes at STR loci.

    First yield is a tuple of names of the fields in details.
    Every subsequent yield is described in the yields section below.

    Parameters
    ----------
    imputation_run_name:
        which imputation run to load genotypes from?
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
    
        If hardcalls, then instead of that just an array nx2 of length alleles
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

        None if hardcalls.
    locus_details:
        tuple of strings with the same length as the first yield
        with the corresponding order.

        None if hardcalls.
        

    Notes
    -----
    Hardcalls mentioned in the locus details are phased hardcalls, and in some
    corner cases will not correspond to the maximum likelihood unphased allele.
    """

    if hardcalls:
        assert var_subset is not None and not details

    chrom, region_poses = region.split(':')
    region_start, _ = [int(pos) for pos in region_poses.split('-')]
    vcf_fname = (f'{ukb}/str_imputed/runs/{imputation_run_name}/'
                f'vcfs/annotated_strs/chr{chrom}.vcf.gz')
    vcf = cyvcf2.VCF(vcf_fname)

    if details:
        yield (
            'motif',
            'period',
            'ref_len',
            'total_per_allele_dosages',
            'total_hardcall_alleles',
            'total_hardcall_genotypes',
            'subset_total_per_allele_dosages',
            'subset_total_hardcall_alleles',
            'subset_total_hardcall_genotypes',
            'subset_het',
            'subset_entropy',
            'subset_HWEP',
            'subset_allele_dosage_r2'
        )
    for record in vcf(region):
        if record.POS < region_start:
            # records that overlap this region but started before this region
            # should be considered part of the pervious region and not returned here
            continue
        if record.INFO.get('PERIOD') is None:
            # there are a few duplicate loci which I didn't handle
            # properly, this identifies and removes them
            continue

        if var_subset is not None and record.POS not in var_subset:
            continue

        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype='beagle-hipstr')

        len_alleles = [trrecord.ref_allele_length] + trrecord.alt_allele_lengths
        len_alleles = [round(allele_len, allele_len_precision) for allele_len in len_alleles]

        if hardcalls:
            yield (trrecord.GetLengthGenotypes()[samples, :-1], np.unique(len_alleles), trrecord.chrom, trrecord.pos, None, None)
            continue

        if details:
            total_dosages = {_len: 0 for _len in np.unique(len_alleles)}
            for p in (1, 2):
                ap = trrecord.format[f'AP{p}']
                total_dosages[len_alleles[0]] += np.sum(np.maximum(0, 1 - np.sum(ap, axis=1)))
                for i in range(ap.shape[1]):
                    total_dosages[len_alleles[i+1]] += np.sum(ap[:, i])

            total_hardcall_alleles = clean_len_alleles(trrecord.GetAlleleCounts())
            total_hardcall_genotypes = clean_len_allele_pairs(trrecord.GetGenotypeCounts())

        if isinstance(samples, slice):
            assert samples == slice(None)
            n_subset_samples = trrecord.GetNumSamples()
        else:
            n_subset_samples = int(np.sum(samples))

        subset_dosage_gts = {
            _len: np.zeros((n_subset_samples, 2)) for _len in np.unique(len_alleles)
        }

        for p in (1, 2):
            # todo genotype dosages
            ap = trrecord.format[f'AP{p}']
            subset_dosage_gts[len_alleles[0]][:, (p-1)] += \
                    np.maximum(0, 1 - np.sum(ap[samples, :], axis=1))
            for i in range(ap.shape[1]):
                subset_dosage_gts[len_alleles[i+1]][:, (p-1)] += ap[samples, i]

        subset_total_dosages = {
            _len: np.sum(subset_dosage_gts[_len]) for _len in subset_dosage_gts
        }

        if details:
            subset_total_hardcall_alleles = clean_len_alleles(trrecord.GetAlleleCounts(samples))
            subset_total_hardcall_genotypes = clean_len_allele_pairs(trrecord.GetGenotypeCounts(samples))
            subset_hardcall_allele_freqs = clean_len_alleles(trrecord.GetAlleleFreqs(samples))

            subset_het = utils.GetHeterozygosity(subset_hardcall_allele_freqs)
            subset_entropy = utils.GetEntropy(subset_hardcall_allele_freqs)
            subset_hwep = utils.GetHardyWeinbergBinomialTest(
                subset_hardcall_allele_freqs,
                subset_total_hardcall_genotypes
            )

            # https://www.cell.com/ajhg/fulltext/S0002-9297(09)00012-3#app1
            # Browning, Brian L., and Sharon R. Browning. "A unified approach to genotype imputation and haplotype-phase inference for large data sets of trios and unrelated individuals." The American Journal of Human Genetics 84.2 (2009): 210-223.
            # appendix 1
            subset_allele_dosage_r2 = {}

            subset_hardcalls = np.around(trrecord.GetLengthGenotypes()[samples, :-1], allele_len_precision)
            for length in len_alleles:
                # calculate allele dosage r**2 for this length
                if length in subset_allele_dosage_r2:
                    continue

                calls = subset_hardcalls == length

                subset_allele_dosage_r2[length] = np.corrcoef(
                    calls.reshape(-1), subset_dosage_gts[length].reshape(-1)
                )[0,1]**2

            locus_details = (
                trrecord.motif,
                str(len(trrecord.motif)),
                str(round(trrecord.ref_allele_length, allele_len_precision)),
                dict_str(round_vals(total_dosages, dosage_precision)),
                dict_str(total_hardcall_alleles),
                dict_str(total_hardcall_genotypes),
                dict_str(round_vals(subset_total_dosages, dosage_precision)),
                dict_str(subset_total_hardcall_alleles),
                dict_str(subset_total_hardcall_genotypes),
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
                trrecord.pos,
                'MAC<20',
                locus_details
            )
            continue

        yield (
            subset_dosage_gts,
            np.unique(len_alleles),
            trrecord.chrom,
            trrecord.pos,
            None,
            locus_details
        )

def filtered_microarray_snps(region):
    """
    Iterate over a region returning genotypes at SNP loci.

    Returns microarray measured SNPs that have been QCed and phased
    (i.e. the ones in the bgen files) which does not include
    all microarray measured SNPs

    Yields tuples
	(genotypes, chrom, pos, alleles, locus_filtered)
    currently, locus_filtered is always None

    genotypes are pairs of indicators:
        0: 1st allele
        1: 2nd allele
    note that bgen files don't distinguish between reference and alt
    alleles. 1st allele may or may not be the reference. (This
    is different than imputed SNPs!)

    Expecting all samples in the input file to have 0 or 1 probs
    so not filtering any calls

    No filtering for samples with rare genotypes
    Loci with only single genotypes remaining are returned - it is up to
    calling code to filter these out
    """
    chrom, posses = region.split(':')
    start, end = posses.split('-')
    start = int(start)
    end = int(end)
    bgen_fname = f'{ukb}/microarray/ukb_hap_chr{chrom}_v2.bgen'
    bgen = bgen_reader.open_bgen(bgen_fname,
                                 allow_complex=True,
                                 verbose=False)

    for variant_num, pos in enumerate(bgen.positions):
        if pos < start:
            continue
        if pos > end:
            break

        probs, missing, ploidy = bgen.read(variant_num,
                                           return_missings=True,
                                           return_ploidies=True)

        # make sure the record looks as expected
        assert np.all(np.logical_or(probs == 0, probs == 1))
        assert probs.shape[2] == 4
        assert bgen.phased[variant_num]
        assert np.all(ploidy == 2)
        assert not np.any(missing)

        # Since all alleles are bialleleic, the genotype can
        # be written as the presence or absence of the second allele
        # for each haplotype for each participant
        reformatted_gts = probs[:, 0, [1,3]]

        yield (
            reformatted_gts,
            chrom,
            pos,
            bgen.allele_ids[variant_num],
            None,
            None
        )


def load_imputed_snps(region: str,
                      samples: np.ndarray,
                      info_thresh: Optional[float] = None,
                      apply_filter: bool = True,
                      details: bool = True,
                      var_subset: Optional[Set[Tuple[int, str, str]]] = None,
                      hardcalls=False): # pos, ref, alt
    """
    Iterate over a region returning genotypes at SNP loci.

    First yield is a tuple of names of the fields in details.
    Every subsequent yield is described in the yields section below.

    Parameters
    ----------
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

        If hardcalls, then instead of that just an array nx2 of length alleles
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

        None if hardcalls
    locus_details:
        tuple of strings with the same length as the first yield
        with the corresponding order.

        None if hardcalls

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

    if hardcalls:
        assert var_subset is not None and not details

    chrom, posses = region.split(':')
    start, end = tuple(int(val) for val in posses.split('-'))
    bgen_fname = f'{ukb}/array_imputed/ukb_imp_chr{chrom}_v3.bgen'
    mfi_fname = f'{ukb}/array_imputed/ukb_mfi_chr{chrom}_v3.txt'
    bgen = bgen_reader.open_bgen(bgen_fname,
                                 verbose=False,
                                 allow_complex=True)

    if details:
        yield (
            'info',
            "total_per_allele_dosages",
            'total_hardcalls',
            'subset_total_per_allele_dosages',
            'subset_total_hardcalls',
            'subset_HWEP'
        )

    with open(mfi_fname) as mfi:
        for variant_num, pos in enumerate(bgen.positions):
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

            if hardcalls:
                probs, missing, ploidy = bgen.read(
                    variant_num,
                    return_missings=True,
                    return_ploidies=True
                )
                probs = np.squeeze(probs)[samples, :]
                gts = np.argmax(probs, axis=1)
                yield (gts, alleles, chrom, pos, None, None)
                continue

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
            assert len(probs.shape) == 2
            assert probs.shape[1] == 3
            assert not bgen.phased[variant_num]
            assert np.all(ploidy == 2)
            assert not np.any(missing)

            dosages = probs[:, 1] + 2*probs[:, 2]
            if details:
                total_alt_dosage = np.sum(dosages)
                total_ref_dosage = 2*dosages.shape[0] - total_alt_dosage

            subset_dosages = dosages[samples]
            subset_total_alt_dosage = np.sum(subset_dosages)
            subset_total_ref_dosage = \
                    2*subset_dosages.shape[0] - subset_total_alt_dosage

            if details:
                hardcalls = probs.argmax(axis=1)
                n_hom_ref = np.sum(hardcalls == 0)
                n_het = np.sum(hardcalls == 1)
                n_hom_alt = np.sum(hardcalls == 2)

                subset_hardcalls = hardcalls[samples]
                subset_n_hom_ref = np.sum(subset_hardcalls == 0)
                subset_n_het = np.sum(subset_hardcalls == 1)
                subset_n_hom_alt = np.sum(subset_hardcalls == 2)

                subset_n_ref_alleles = 2*subset_n_hom_ref + subset_n_het
                subset_frac_ref_alleles = (subset_n_ref_alleles)/(2*subset_hardcalls.shape[0])
                subset_exp_hom_frac = subset_frac_ref_alleles**2 + (1-subset_frac_ref_alleles)**2
                subset_hwep = scipy.stats.binom_test(subset_n_hom_ref + subset_n_hom_alt,
                                                     n=subset_hardcalls.shape[0],
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

            if (subset_total_alt_dosage < 20 or subset_total_ref_dosage < 20) and apply_filter:
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

