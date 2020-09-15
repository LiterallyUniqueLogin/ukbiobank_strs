#!/bin/env python3

"""
Read VCFs of STRs or SNPs,
apply filters to them for quality control
and then return yield them
"""

import os
import os.path

import cyvcf2
import numpy as np

ukb = os.environ['UKB']


def filtered_strs(imputation_run_name,
                       region):
    """
    Iterate over a region returning genotypes at STR loci.

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

