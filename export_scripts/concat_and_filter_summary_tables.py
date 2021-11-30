#!/usr/bin/env python3

import argparse
import os

import numpy as np

ukb = os.environ['UKB']

def write_line(out, split, indices, extra):
    out.write(
        extra + '\t' + '\t'.join(split[idx] for idx in indices) + '\n'
    )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotypes', nargs='+')
    args = parser.parse_args()

    phenotypes = args.phenotypes

    col_descs = (
        'phenotype\n'
        'chrom\n'
        'start_pos - the start bp of the STR, 1-based\n'
        'repeat_unit - as preliminarily inferred by TRTools from the reference STR '
        'sequence and the repeat period stated in the SNPSTR reference panel\n'
        'multiallelicness - amongst the population being tested for association, '
        'this is the fraction of total allelic dosage at this locus '
        'taken up by all alleles but the two most common alleles\n'
        'association_p_value\n'
        "pcausal - FINEMAP's posterior probability of causality\n"
        'relation_to_gene - if this STR is transcribed, for each transcript '
        'what is the GENCODE gene type of that transcript (i.e. protein coding, '
        'lncRNA, etc.) and what is the GENCODE feature type of the region the STR '
        'is in (i.e. intron, exon, etc.)\n'
        "transcribed - each transcript this STR is in (if any), including the transcript's "
        "GENCODE gene type (i.e. protein coding, lncRNA, etc.) and the transcript's "
        'transcript support level (1-5 or missing)\n'
    )

    with open(f'{ukb}/export_scripts/results/putatively_causal_STRs_README.txt', 'w') as readme:
        readme.write(
            'Any STR x phenotype association that has p-value <= 1e-10 and '
            'FINEMAP posterior probability of causality >= 0.8\n'
            + col_descs
        )

    with open(f'{ukb}/export_scripts/results/curated_STRs_README.txt', 'w') as readme:
        readme.write(
            'Any curated STR x phenotype association\n'
            + col_descs
        )

    cols = [
        'chrom',
        'start_pos',
        'repeat_unit',
        'subset_multiallelicness',
        'association_p_value',
        'pcausal',
        'relation_to_gene',
        'transcribed'
    ]

    with open(f'{ukb}/export_scripts/results/putatively_causal_STRs.tab', 'w') as strong_out, \
            open(f'{ukb}/export_scripts/results/curated_STRs.tab', 'w') as curated_out:
        first_phen = True
        for phenotype in phenotypes:
            with open(f'{ukb}/finemapping/summary/{phenotype}_table.tab') as table:
                if first_phen:
                    header = next(table)
                    split = header.strip().split('\t')
                    indices = [split.index(col) for col in cols]
                    association_p_value_index = split.index('association_p_value')
                    pcausal_index = split.index('pcausal')
                    subset_multiallelicness_index = split.index('subset_multiallelicness')
                    curated_index = split.index('curated')
                    only_curated_index = split.index('included_only_due_to_curation')
                    only_lit_index = split.index('included_only_due_to_literature')
                    split[subset_multiallelicness_index] = 'multiallelicness'
                    write_line(strong_out, split, indices, 'phenotype')
                    write_line(curated_out, split, indices, 'phenotype')
                else:
                    next_line = next(table)
                    if next_line != header:
                        print(next_line, header)
                        assert False

                for line in table:
                    split = line.strip().split('\t')
                    if split[curated_index] == 'True':
                        write_line(curated_out, split, indices, phenotype)
                    association_p_value = float(split[association_p_value_index])
                    if association_p_value != 0 and -np.log10(association_p_value) < 10:
                        continue
                    if split[pcausal_index] != 'NA' and float(split[pcausal_index]) < .8:
                        continue
                    if split[only_curated_index] == 'True' or split[only_lit_index] == 'True':
                        continue
                    write_line(strong_out, split, indices, phenotype)
            first_phen = False

if __name__ == '__main__':
    main()
