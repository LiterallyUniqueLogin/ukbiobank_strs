#!/bin/env python3

"""
Read VCFs of STRs or SNPs,
apply filters to them for quality control
and then return yield them
"""

import os
import os.path
import time
from typing import List, Tuple
import sys

import bgen_reader
import cyvcf2
import numpy as np
import numpy.ma

ukb = os.environ['UKB']


def filtered_strs(imputation_run_name,
                  region):
    """
    Iterate over a region of STRs returning genotypes at each locus.

	First value yielded is the array of samples in the VCF
	Every subsequent value is the tuple
	(genotypes, chrom, pos, alleles, locus_info, locus_filtered)
    currently, locus_filtered is always None

    genotypes are pairs of haplotypes
    measured by difference in # repeats from the reference

    Samples with imputed genotypes with too low an expected probability
    (as measured by product of their allele probabilities)
    have their genotypes set to nan

    No filtering for samples with rare genotypes
    Loci with only single genotypes remaining are returned - it is up to
    calling code to filter these out
    """
    chrom, _ = region.split(':')
    vcf_fname = (f'{ukb}/str_imputed/runs/{imputation_run_name}/'
                f'vcfs/annotated_strs/chr{chrom}.vcf.gz')
    vcf = cyvcf2.VCF(vcf_fname)
    yield np.array(vcf.samples)

    vcf_region = vcf(region)
    for locus in vcf_region:
        if 'PERIOD' not in dict(locus.INFO):
            # there are a few duplicate loci which I didn't handle
            # properly, this identifies and removes them
            continue
        # setup for this locus

        # gt entries are sequence allele indexes
        idx_gts = locus.genotype.array()[:, :2]
        seq_allele_lens = [len(allele) for allele in locus.ALT]
        seq_allele_lens.insert(0, len(locus.REF))

        # modify gt entries to be length alleles
        len_gts = np.zeros(idx_gts.shape, dtype=float)
        for seq_allele_idx, seq_allele_len in enumerate(seq_allele_lens):
            len_gts[idx_gts == seq_allele_idx] = seq_allele_len

        # get Beagle genotype probability
        aps = []
        for copy in 1, 2:
            ap = locus.format(f'AP{copy}')
            ref_prob = np.maximum(0, 1 - np.sum(ap, axis=1))
            ref_prob = ref_prob.reshape(-1, 1)
            ap = np.concatenate((ref_prob, ap), axis=1)
            ap = ap[np.arange(ap.shape[0]), idx_gts[:, (copy - 1)]]
            aps.append(ap)
        gp = np.multiply(aps[0], aps[1])

        # filter calls whose phased genotype probability
        # as estimated by Beagle is below .9
        # convert to float so we can use np.nan
        filtered_samples = gp < .9
        len_gts[filtered_samples, :] = np.nan
        n_filtered_samples = np.sum(filtered_samples)

        # convert from length in bp to length in diff from ref in repeat
        # units
        len_gts -= seq_allele_lens[0]
        len_gts /= locus.INFO['PERIOD']

        yield (
            len_gts,
            locus.CHROM,
            locus.POS,
            ','.join(str(l) for l in np.unique(len_gts[~np.isnan(len_gts)])),
            f'PERIOD={locus.INFO["PERIOD"]},REF={locus.REF}',
            None
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


def filtered_imputed_snps(region):
    """
    Iterate over a region returning genotypes at SNP loci.

    Returns imputed SNPs

    Yields tuples
	(genotypes, chrom, pos, allele, locus_filtered)
    Exactly one of genotypes and locus_filtered will be None

    genotypes are pairs of indicators:
        0: 1st allele
        1: 2nd allele
    1st allele is always the reference for imputed SNPs.
    (This is different than microarray SNPs!)
    Also note that these genotypes are not phased so the ordering
    of the genotypes for heterozygous samples is arbitrary.

    Loci with INFO < 0.3 or so should be filtered
    But currently are not to allow this to be tuned downstream.
    Samples where the probability of their most probable genotypes is <.9
    are filtered

    No filtering for samples with rare genotypes
    Loci with only single genotypes remaining are returned - it is up to
    calling code to filter these out
    """
    chrom, posses = region.split(':')
    start, end = posses.split('-')
    start = int(start)
    end = int(end)
    bgen_fname = f'{ukb}/array_imputed/ukb_imp_chr{chrom}_v3.bgen'
    mfi_fname = f'{ukb}/array_imputed/ukb_mfi_chr{chrom}_v3.txt'
    bgen = bgen_reader.open_bgen(bgen_fname,
                                 verbose=False,
                                 allow_complex=True)

    with open(mfi_fname) as mfi:
        for variant_num, pos in enumerate(bgen.positions):
            mfi_line = next(mfi)
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
                    'info=NA',
                    'MAF=0'
                )
            '''
            elif float(info_str) < 0.3:
                yield (
                    None,
                    chrom,
                    pos,
                    bgen.allele_ids[variant_num],
                    f'info={info_str}',
                    'info<0.3'
                )
                continue
            '''

            probs, missing, ploidy = bgen.read(variant_num,
                                           return_missings=True,
                                           return_ploidies=True)
            probs = np.squeeze(probs)

            # make sure the record looks as expected
            assert probs.shape[1] == 3
            assert not bgen.phased[variant_num]
            assert np.all(ploidy == 2)
            assert not np.any(missing)

            best_genotype = probs.argmax(axis=1)
            highest_prob = probs.max(axis=1)
            out_gts = np.zeros((probs.shape[0], 2))

            out_gts[best_genotype == 2, 0]  = 1
            out_gts[best_genotype > 0, 1]  = 1
            out_gts[highest_prob < .9, :] = np.nan


            yield (
                out_gts,
                chrom,
                pos,
                bgen.allele_ids[variant_num],
                f'info={info_str}',
                None
            )


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


