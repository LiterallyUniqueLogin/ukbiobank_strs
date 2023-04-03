#!/usr/bin/env python3

import argparse
import shutil

import cyvcf2
import numpy as np
import scipy.stats

import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

def main():
    parser = argparse.ArgumentParser()
    # f'{ukb}/export_scripts/intermediate_results/chr{chrom}_loci_summary.tab'
    parser.add_argument('out')
    # vcf_fname = (f'{ukb}/str_imputed/runs/{str_imputation_run_name}/'
    #                         f'vcfs/annotated_strs/chr{chrom}.vcf.gz')
    parser.add_argument('vcf_fname')
    parser.add_argument('all_white_brits_fname')
    args = parser.parse_args()

    with open(f'{args.out}.temp', 'w') as output:
        main_helper(output, args.vcf_fname, args.all_white_brits_fname)

    shutil.move(f'{args.out}.temp', args.out)

def main_helper(output, vcf_fname, all_white_brits_fname):
    vcf = cyvcf2.VCF(vcf_fname)

    subset_samples = []
    with open(all_white_brits_fname) as samp_file:
        next(samp_file)
        for line in samp_file:
            subset_samples.append(line.strip())

    all_samples = [l[0] for l in np.char.split(vcf.samples, '_')]
    samples = np.isin(all_samples, subset_samples)

    output.write('chr\tpos\tallele_dist\tentropy\theterozygosity\tmultiallelicness\n')
    for record in vcf:
        if record.INFO.get('PERIOD') is None:
            continue
        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype='hipstr')

        len_alleles = [trrecord.ref_allele_length] + trrecord.alt_allele_lengths
        total_subset_dosages = {
            len_: 0 for len_ in np.unique(len_alleles)
        }
        for p in (1, 2):
            ap = trrecord.format[f'AP{p}']
            total_subset_dosages[len_alleles[0]] += \
                    np.sum(np.maximum(0, 1 - np.sum(ap[samples, :], axis=1)))
            for i in range(ap.shape[1]):
                total_subset_dosages[len_alleles[i+1]] += np.sum(ap[samples, i])

        for len_ in total_subset_dosages:
            total_subset_dosages[len_] /= np.sum(samples)*2

        entropy = scipy.stats.entropy(list(total_subset_dosages.values()), base=2)
        heterozygosity = 1 - np.sum(val**2 for val in total_subset_dosages.values())
        multiallelicness = sum(sorted(total_subset_dosages.values())[:-2])

        output.write(
            trrecord.chrom + '\t' + str(trrecord.pos) + '\t' +
            str(total_subset_dosages) + '\t' + str(entropy) + '\t' +
            str(heterozygosity) + '\t' + str(multiallelicness) + '\n'
        )

if __name__ == '__main__':
    main()
