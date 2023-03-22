#!/usr/bin/env python3

import argparse

import cyvcf2
import numpy as np

def check(merged_file, batches):
    merged_samples = np.array(merged_file.samples)
    batch_samples = [np.array(batch.samples) for batch in batches]

#    if not np.all(merged_samples[:-1] < merged_samples[1:]):
#        print(merged_samples)
#        assert False
#
#    merged_reformatted_samples = np.char.add(np.char.add(merged_samples, '_'), merged_samples)
#    for batch_index, batch_sample_list in enumerate(batch_samples):
#        if not np.all(batch_sample_list[:-1] < batch_sample_list[1:]):
#            print(batch_sample_list)
#            assert False

    if not np.all(merged_samples == np.concatenate(batch_samples)):
        print(repr(merged_samples))
        print(repr(batch_samples))
        assert False

    for merged_var, batch_vars in zip(merged_file, zip(*batches)):
        merged_gts = merged_var.genotype.array()[:, :-1]
        merged_alleles = [merged_var.REF] + merged_var.ALT
        sample_offset = 0
        for batch_idx, batch_var in enumerate(batch_vars):
            if merged_var.INFO['START'] != batch_var.INFO['START'] or merged_var.INFO['END'] != batch_var.INFO['END']:
                print(merged_var.POS, batch_idx)
                assert False

            batch_gts = batch_var.genotype.array()[:, :-1]
            corresponding_merged_gts = merged_gts[sample_offset:(sample_offset+len(batch_samples[batch_idx])), :]

            len_offset = len(batch_var.REF) - batch_var.INFO['END'] + batch_var.INFO['START'] - 1
            batch_alleles = [batch_var.REF] + batch_var.ALT
            for allele_idx, allele in enumerate(batch_alleles):
                if np.sum(batch_gts == allele_idx) == 0:
                    continue
                merged_allele_set = np.unique(corresponding_merged_gts[batch_gts == allele_idx])
                if not len(merged_allele_set) == 1:
                    print(merged_var.POS, batch_idx, allele_idx)
                    assert False
                merged_allele = merged_alleles[merged_allele_set[0]]
                if not len(allele) - len_offset == len(merged_allele):
                    print(merged_var.POS, batch_idx, allele_idx)
                    assert False

            sample_offset += len(batch_samples[batch_idx])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('merged_file')
    parser.add_argument('batches', nargs='+')
    args = parser.parse_args()

    try:
        batch_files = []
        print("opening merged file")
        merged_file = cyvcf2.VCF(args.merged_file)
        for i, file in enumerate(args.batches):
            print(f"opening batch file {i+1}")
            batch_files.append(cyvcf2.VCF(file))

        check(merged_file, batch_files)
    finally:
        merged_file.close()
        for file in batch_files:
            file.close()

    print('Success')

if __name__ == "__main__":
    main()
