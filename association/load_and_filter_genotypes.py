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
    Iterate over a region returning genotypes at STR loci.

	First value yielded is the array of samples in the VCF
	Every subsequent value is the tuple
	(genotypes, chrom, pos, filtered_samples, filtered_rare_alleles)

    genotypes have the unit: difference in # repeats from the reference

    Samples with imputed genotypes with too low an expected probability
    have their genotypes set to nan
    Samples with genotypes that are too rare have their genotypes set to nan

    If after filtering a locus has only one genotype, it is skipped
    and then next qualifying locus is yielded
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

        # filter length alleles with too few occurrences
        len_alleles, counts = np.unique(
            len_gts[~np.isnan(len_gts)],
            return_counts=True
        )
        filtered_rare_alleles = 0
        for len_allele_idx, len_allele in enumerate(len_alleles):
            if (counts[len_allele_idx] > 0 and
                    counts[len_allele_idx] < 5):
                len_gts[len_gts == len_allele] = np.nan
                filtered_rare_alleles += 1

        # make sure this locus doesn't have only one nonfiltered genotype
        # (or worse, no nonfiltered genotypes)
        not_nans = len_gts[~np.isnan(len_gts)]
        trivial_data = len(not_nans) == 0 or np.all(not_nans == not_nans[0])
        if trivial_data:
            continue

        # convert from length in bp to length in diff from ref in repeat
        # units
        len_gts -= seq_allele_lens[0]
        len_gts /= locus.INFO['PERIOD']

        yield (
            len_gts,
            locus.CHROM,
            locus.POS,
            n_filtered_samples,
            filtered_rare_alleles
        )


def filtered_snps(region):
    """
    Iterate over a region returning genotypes at SNP loci.

    Currently returns microarray measured SNPs that have been QCed,
    not all microarray measured SNPs, and not imputed SNPs

	First value yielded is the array of samples in the VCF
	Every subsequent value is the tuple
	(genotypes, chrom, pos, filtered_samples, filtered_rare_alleles)

    genotypes have the unit:
        0 (hom. 1st allele),
        0.5 (het)
        1 (hom. 2nd allele)
    note that bgen files don't distinguish between reference and alt
    alleles. 1st allele may or may not be the reference

    Expecting all samples in the input file to have 0 or 1 probs
    so filtered_samples will always be zero

    If either allele is too rare, the locus will be skipped for being
    too homozygous

    Because of this, any locus that is yielded will ahve filtered_rare_alleles
    as zero
    """
    chrom, posses = region.split(':')
    start, end = posses.split('-')
    start = int(start)
    end = int(end)
    bgen_fname = f'{ukb}/microarray/ukb_hap_chr{chrom}_v2.bgen'
    samples_fname = f'{ukb}/microarray/ukb46122_hap_chr1_v2_s487314.sample'
    bgen = bgen_reader.read_bgen(bgen_fname,
                                 samples_filepath=samples_fname,
                                 verbose=False)
    variants = bgen["variants"].compute()
    genotypes = bgen["genotype"]
    samples = bgen["samples"].to_numpy(dtype='U7')
    samples = np.char.add(np.char.add(samples, '_'), samples)
    yield samples

    for variant_num, pos in enumerate(variants['pos']):
        if pos < start:
            continue
        if pos > end:
            break

        genotype_info = genotypes[variant_num].compute()

        # make sure the record looks as expected
        assert np.all(np.logical_or(
            genotype_info['probs'] == 0, genotype_info['probs'] == 1
        ))
        assert genotype_info['probs'].shape[1] == 4
        assert genotype_info['phased']
        assert np.all(genotype_info['ploidy'] == 2)
        assert not np.any(genotype_info['missing'])

        # for this dataset, these aren't probs but hardcalls
        gts = genotype_info['probs']

        n_first_allele = np.sum(gts[:, [0,2]])
        if n_first_allele <  5:
            continue
        n_second_allele = np.sum(gts[:, [1,3]])
        if n_second_allele <  5:
            continue

        # Since all alleles are bialleleic, the genotype can
        # be written as the presence or absence of the second allele
        # for each haplotype for each participant
        reformatted_gts = gts[:, [1,3]]

        yield (
            reformatted_gts,
            chrom,
            pos,
            0,
            0
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


