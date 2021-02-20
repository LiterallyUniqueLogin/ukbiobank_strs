#!/bin/env python3

"""
Read VCFs of STRs or SNPs,
apply filters to them for quality control
and then return yield them
"""

import os
import os.path
import time
from typing import List, Optional, Tuple
import sys

import bgen_reader
import cyvcf2
import numpy as np
import numpy.ma
import scipy

ukb = os.environ['UKB']

sys.path.insert(0, f'{ukb}/../trtools/repo')

import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

def load_strs(imputation_run_name: str,
              region: str,
              samples: np.ndarray):
    """
    Iterate over a region returning genotypes at SNP loci.

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
    gt:
        A 1D array of length n-samples
        Each value is the length dosage measured in number of repeats
        different from the reference.

        None if locus_filtered is not None
    chrom: str
        e.g. '13'
    pos: int
    alleles: str
        e.g. 'ACAC,ACACAC'
    locus_filtered:
        None if the locus is not filtered, otherwise
        a string explaining why.
        'MAC<20' if there are fewer than 20 minor allele hardcalls
        after sample subsetting, per plink's standard
    locus_details:
        tuple of strings with the same length as the first yield
        with the corresponding order.
    """

    chrom, _ = region.split(':')
    vcf_fname = (f'{ukb}/str_imputed/runs/{imputation_run_name}/'
                f'vcfs/annotated_strs/chr{chrom}.vcf.gz')
    vcf = cyvcf2.VCF(vcf_fname)
   
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
        if record.INFO.get('PERIOD') is None:
            # there are a few duplicate loci which I didn't handle
            # properly, this identifies and removes them
            continue

        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype='beagle-hipstr')

        len_alleles = [trrecord.ref_allele_length] + trrecord.alt_allele_lengths

        total_dosages = {_len: 0 for _len in np.unique(len_alleles)}
        for p in (1, 2):
            ap = trrecord.format[f'AP{p}']
            total_dosages[len_alleles[0]] += np.sum(np.maximum(0, 1 - np.sum(ap, axis=1)))
            for i in range(ap.shape[1]):
                total_dosages[len_alleles[i+1]] = np.sum(ap[:, i])

        total_hardcall_alleles = trrecord.GetAlleleCounts()
        total_hardcall_genotypes = trrecord.GetGenotypeCounts()

        subset_total_dosages = {_len: 0 for _len in np.unique(len_alleles)}
        for p in (1, 2):
            ap = trrecord.format[f'AP{p}']
            subset_total_dosages[len_alleles[0]] += np.sum(np.maximum(0, 1 - np.sum(ap[samples, :], axis=1)))
            for i in range(ap.shape[1]):
                subset_total_dosages[len_alleles[i+1]] = np.sum(ap[samples, i])

        subset_total_hardcall_alleles = trrecord.GetAlleleCounts(samples)
        subset_total_hardcall_genotypes = trrecord.GetGenotypeCounts(samples)
        subset_hardcall_allele_freqs = trrecord.GetAlleleFreqs(samples)

        subset_het = utils.GetHeterozygosity(subset_hardcall_allele_freqs)
        subset_entropy = utils.GetEntropy(subset_hardcall_allele_freqs)
        subset_hwep = utils.GetHardyWeinbergBinomialTest(
            subset_hardcall_allele_freqs,
            subset_total_hardcall_genotypes
        )

        subset_dosages = trrecord.GetDosages()[samples]

        # https://www.cell.com/ajhg/fulltext/S0002-9297(09)00012-3#app1
        # Browning, Brian L., and Sharon R. Browning. "A unified approach to genotype imputation and haplotype-phase inference for large data sets of trios and unrelated individuals." The American Journal of Human Genetics 84.2 (2009): 210-223.
        # appendix 1
        subset_allele_dosage_r2 = {}

        allele_lens = [trrecord.ref_allele_length] + trrecord.alt_allele_lengths
        subset_hardcall_indicies = trrecord.GetGenotypeIndicies()[samples, :-1]
        aps = []
        for p in 1, 2:
            ap = trrecord.format[f'AP{p}'][samples, :]
            ref_ap = np.maximum(0, np.sum(ap, axis=1))
            aps.append(np.concatenate((ref_ap.reshape(-1, 1), ap), axis=1))
        stacked_aps = np.stack(aps, axis=1)
        assert stacked_aps.shape[0] == subset_dosages.shape[0]
        assert stacked_aps.shape[1] == 2
        for length in allele_lens:
            # calculate allele dosage r**2 for this length
            for other_length in subset_allele_dosage_r2:
                if np.isclose(length, other_length):
                    continue

            calls = np.zeros(subset_hardcall_indicies.shape, dtype=bool)
            for idx, inner_length in enumerate(allele_lens):
                if not np.isclose(inner_length, length):
                    continue
                calls = calls | (subset_hardcall_indicies == idx)

            same_len_alleles = [np.isclose(allele_len, length) for allele_len in allele_lens]
            probs = np.sum(stacked_aps[:, :, same_len_alleles], axis=2)

            assert probs.shape == calls.shape

            subset_allele_dosage_r2[length] = \
                    np.corrcoef(calls.reshape(-1), probs.reshape(-1))[0,1]**2

        locus_details = (
            trrecord.motif,
            str(len(trrecord.motif)),
            str(trrecord.ref_allele_length),
            str(total_dosages),
            str(total_hardcall_alleles),
            str(total_hardcall_genotypes),
            str(subset_total_dosages),
            str(subset_total_hardcall_alleles),
            str(subset_total_hardcall_genotypes),
            str(subset_het),
            str(subset_entropy),
            str(subset_hwep),
            str(subset_allele_dosage_r2)
        )

        mac = list(subset_total_hardcall_alleles.values())
        mac.pop(np.argmax(mac))
        mac_lt_20 = np.sum(mac) < 20

        if mac_lt_20:
            yield (
                None,
                trrecord.chrom,
                trrecord.pos,
                ','.join(str(_len) for _len in np.unique(len_alleles)),
                'MAC<20',
                locus_details
            )
            continue

        yield (
            subset_dosages,
            trrecord.chrom,
            trrecord.pos,
            ','.join(str(_len) for _len in np.unique(len_alleles)),
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
                      return_dosages: bool = True,
                      info_thresh: Optional[float] = None):
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
    return_dosages:
        whether or not to return dosages or hardcalls
    info_thresh:
        Loci must have INFO scores of at least this thresh or will be filtered

    Yields
    ------
    gt:
        Either dosages or hardcalls. If dosages:
        A 1D array of length n-samples
        Each value is the dosage of the alternate allele:
            prob(AB) + 2*prob(BB)

        If hardcalls:
        A 1D array of length n-samples
        0,1,2, corresponding to the whichever of AA,AB or BB
        has max likelihood.
        Entries filtered by call_thresh are np.nan

        None if locus_filtered is not None
    chrom: str
        e.g. '13'
    pos: int
    alleles: str
        e.g. 'A,G'
    locus_filtered:
        None if the locus is not filtered, otherwise
        a string explaining why.
        'MAF=0' if there is only one allele present
        in the calls at this locus.
        'MAC<20' if there are fewer than 20 minor allele hardcalls
        after sample subsetting, per plink's standard
        'info<{thresh}' if an info thresh
        is set and this locus is under that
    locus_details:
        tuple of strings with the same length as the first yield
        with the corresponding order.

    Notes
    -----
    1st allele is always the reference for imputed SNPs, alt allele
    is always nonreference.
    (This is different than microarray SNPs!)

    (Also note that these genotypes are not phased so the ordering
    of the genotypes for heterozygous samples is arbitrary. This
    currently doesn't matter as hard calls are)
    """
    chrom, posses = region.split(':')
    start, end = tuple(int(val) for val in posses.split('-'))
    bgen_fname = f'{ukb}/array_imputed/ukb_imp_chr{chrom}_v3.bgen'
    mfi_fname = f'{ukb}/array_imputed/ukb_mfi_chr{chrom}_v3.txt'
    bgen = bgen_reader.open_bgen(bgen_fname,
                                 verbose=False,
                                 allow_complex=True)

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

            info_str = mfi_line.split()[-1]
            if info_str == 'NA':
                yield (
                    None,
                    chrom,
                    pos,
                    bgen.allele_ids[variant_num],
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
            total_alt_dosage = np.sum(dosages)
            total_ref_dosage = 2*dosages.shape[0] - total_alt_dosage

            subset_dosages = dosages[samples]
            subset_total_alt_dosage = np.sum(subset_dosages)
            subset_total_ref_dosage = \
                    2*subset_dosages.shape[0] - subset_total_alt_dosage

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

            if (subset_n_ref_alleles < 20 or
                    subset_n_ref_alleles > subset_hardcalls.shape[0]*2 - 20):
                yield (
                    None,
                    chrom,
                    pos,
                    bgen.allele_ids[variant_num],
                    'MAC<20',
                    locus_details
                )
                continue

            if info_thresh is not None and float(info_str) < info_thresh:
                yield (
                    None,
                    chrom,
                    pos,
                    bgen.allele_ids[variant_num],
                    'info<{info_thresh}',
                    locus_details
                )
                continue

            if return_dosages:
                out_gts = subset_dosages
            else:
                out_gts = subset_hardcalls
            yield (
                out_gts,
                chrom,
                pos,
                bgen.allele_ids[variant_num],
                None,
                locus_details
            )


"""
def load_all_haplotypes(variant_generator, samplelist):
    '''
    Load all the haplotypes in an interator.

    Takes in a generator of haplotypes
    which returns 2D arrays of size nsamples x 2
    and returns a 2D array of size (2 x nsamples) x loci

    Parameters
    ----------
    samplelist:
        list of samples to keep

    Returns
    -------
    Tuple[np.ndarray, List[int]]:
        The 2D array of haplotypes and a
        list of positions for each locus
    '''
    gts_per_locus_list = []
    poses = []
    samples = next(variant_generator)

    print("Loading haplotypes ... ", flush=True)
    nloci = 0
    batch_time = 0
    batch_size = 50
    total_time = 0
    start_time = time.time()
    for len_gts, _, pos, _, _  in variant_generator:
        gts_per_locus_list.append(len_gts)
        poses.append(pos)

        nloci += 1
        duration = time.time() - start_time
        total_time += duration
        batch_time += duration
        if nloci % batch_size == 0:
            print(
                f"time/locus (last {batch_size}): "
                f"{batch_time/batch_size:.2e}s\t\t\t"
                f"time/locus ({nloci} total loci): {total_time/nloci:.2e}s",
                flush=True,
                end='\r'
            )
            batch_time = 0
        start_time = time.time()
    sys.stdout.write('\033[2K\033[1G') # erase and go to beginning of line?
    # https://stackoverflow.com/questions/5290994/remove-and-replace-printed-items
    print(
        f"Done loading haplotypes. "
        f"time: {total_time:.2e}s ({nloci} total loci)",
        flush=True
    )
    gts_per_locus = np.stack(gts_per_locus_list, axis=-1)
    samples_to_use = np.isin(samples, samplelist)
    gts_per_locus = gts_per_locus[samples_to_use, :, :]

    # shape is now participant x haplotype x locus
    # all haplotypes are equally useful here, regardless of which person
    # they come from, so combine the first two dimensions into one
    gts_per_locus = gts_per_locus.reshape(gts_per_locus.shape[0]*2,
                                          gts_per_locus.shape[2])
    gts_per_locus = numpy.ma.masked_invalid(gts_per_locus)
    return gts_per_locus, np.array(poses), samples[samples_to_use]

"""
