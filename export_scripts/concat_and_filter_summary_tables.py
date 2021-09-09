#!/usr/bin/env python3

import argparse
import os

import numpy as np

ukb = os.environ['UKB']

def write_line(out, split, indices, extra):
    out.write(
        '\t'.join(split[idx] for idx in indices)
        + '\t' + extra + '\n'
    )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotypes', nargs='+')
    args = parser.parse_args()

    phenotypes = args.phenotypes

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

    with open(f'{ukb}/export/combined/finemapped_loci.tab', 'w') as out:
        first_phen = True
        for phenotype in phenotypes:
            with open(f'{ukb}/export/{phenotype}/{phenotype}_table.tab') as table:
                if first_phen:
                    header = next(table)
                    split = header.strip().split('\t')
                    indices = [split.index(col) for col in cols]
                    association_p_value_index = split.index('association_p_value')
                    pcausal_index = split.index('pcausal')
                    subset_multiallelicness_index = split.index('subset_multiallelicness')
                    write_line(out, split, indices, 'phenotype')
                else:
                    assert next(table) == header

                for line in table:
                    split = line.strip().split('\t')
                    association_p_value = float(split[association_p_value_index])
                    if association_p_value != 0 and -np.log10(association_p_value) < 15:
                        continue
                    if float(split[pcausal_index]) < .4:
                        continue
                    if float(split[subset_multiallelicness_index]) < .05:
                        continue
                    write_line(out, split, indices, phenotype)
            first_phen = False

if __name__ == '__main__':
    main()
