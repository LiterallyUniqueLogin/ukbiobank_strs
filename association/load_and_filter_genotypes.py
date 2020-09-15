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

    Samples with imputed genotypes with too low an expected probability
    have their genotypes set to nan
    Samples with genotypes that are too rare have their genotypes set to nan

    If after filtering a locus has only one genotype, it is skipped
    and then next qualifying locus is yielded
    """
    chrom, _ = region.split(':')
    vcf_floc = (f'{ukb}/str_imputed/runs/{imputation_run_name}/'
                f'vcfs/annotated_strs/chr{chrom}.vcf.gz')
    vcf = cyvcf2.VCF(vcf_floc)
    yield vcf.samples

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
        not_nans = len_gts[~np.isnan(len_gts)]
        trivial_data = np.all(not_nans == not_nans[0])
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

    print("Loading haplotypes ... ")
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
        f"time/locus ({nloci} total loci): {total_time/nloci:.2e}s",
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
    return gts_per_locus, np.array(poses)


