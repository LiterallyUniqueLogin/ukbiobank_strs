import os
import os.path
import sys
import time

import cyvcf2
import numpy as np
import zarr

# Use a LMDB store
# This seems to suggest that LMDB is faster than simply using the filesystem
# https://groups.google.com/forum/#!topic/caffe-users/WqBfwUJ98c0
# Confirm this

# Zarr's default compression algorithm, blosc, uses 8 threads internally. Make
# sure this is respected when submitting this program to the cluster

# create array in Zarr big enough to store all the data we want

zarr_store_loc = 'filtering.lmdb.zarr'
max_n_alleles, n_allele_dtype = 2**8 - 1, 'uint8'


def create_phased_hardcall_store(n_samples, n_variants, chrom):
    with _zarr_store() as store:
        # hardcalls are calls which take the genotype with the highest
        # probability and ignore the probabilities of the rest.
        # phased hardcalls mean that '0|1' and '1|0' are different genotypes
        # and so their probabilities are not combined. These are the GTs that
        # Beagle outputs. Note that it is probably more realistic to work with
        # unphased hardcalls (i.e. p('0/1') which equals p('0|1') + p('1|0'))
        # but calculating those would be more compute intensive and the
        # resulting difference is unlikely to matter

        # Create an array storing the hardcall probabilities
        # Beagle stores AP1 and AP2 as two digit decimals, so AP1*AP2 can be
        # stored as a 4 digit decimal, which fits into two bytes
        # so store these probabilities as 16-bit ints, and to convert
        # from stored value to probability, just divide by 1e4
        # Use (1000,1000) chunks, which equates to 2MB chunks with 16bit ints
        # use -1 as fill value - this is a meaningless probability and so is
        # recognizable as an unfilled location
        zarr.create(
            shape=(n_variants, n_samples),
            chunks=(1000, 1000),
            dtype='int16',
            fill_value=-1,
            store=store,
            path=_phased_hardcall_probs_zpath(chrom)
        )
        # Create an array storing the hardcalls
        # This means that variants can have up to 127 alternate alleles
        # I make no guarantee that the first allele has a lower index than the
        # second allele
        zarr.create(
            shape=(n_variants, n_samples, 2),
            chunks=(1000, 1000, 2),
            dtype=n_allele_dtype,
            fill_value=max_n_alleles - 1,
            store=store,
            path=_phased_hardcalls_zpath(chrom)
        )


def read_phased_hardcalls_from_vcf(vcf_floc, chrom):
    with _zarr_store() as store:
        gts = zarr.open_array(
            store=store,
            path=_phased_hardcalls_zpath(chrom)
        )
        probs = zarr.open_array(
            store=store,
            path=_phased_hardcall_probs_zpath(chrom)
        )
        vcf = cyvcf2.VCF(vcf_floc)

        start = time.time()
        prev = time.time()
        variant_count = 0
        allele_count = 0
        prev_allele_count = 0
        logging_interval = 20  # in variants
        for variant in vcf:
            variant_count += 1
            if variant_count % logging_interval == 0:
                now = time.time()
                interval_time = now - prev
                interval_time_per_variant = interval_time / logging_interval
                interval_time_per_allele = interval_time / prev_allele_count
                total_time = now - start
                total_time_per_variant = total_time / variant_count
                total_time_per_allele = total_time / allele_count
                print(f"Reading variant {variant_count} at POS {variant.POS}")
                print(f"Time per variant last {logging_interval} variants: "
                      f"{interval_time_per_variant:.2f}sec")
                print(f"Time per allele last {logging_interval} "
                      f"variants: {interval_time_per_allele:.2f}sec")
                print(f"Time per variant over all variants: "
                      f"{total_time_per_variant:.2f}sec")
                print(f"Time per allele over all variants: "
                      f"variants: {total_time_per_allele:.2f}sec")
                sys.stdout.flush()
                prev_allele_count = 0
                prev = time.time()

            prev_allele_count += len(variant.ALT) + 1
            allele_count += len(variant.ALT) + 1

            if len(variant.ALT) > max_n_alleles - 1:
                raise OverflowError((
                    f"Variant {variant.POS} had {len(variant.ALT)}alternate "
                    f"alleles but only had space for {max_n_alleles - 1} "
                    f"alternate alleles"
                ))

            # set the genotypes
            gts[variant_count, :, :] = variant.genotype.array()[:, :2]

            # set the probabilities
            best_copy_prob = {}
            for copy in [1, 2]:
                ap = variant.format(f'AP{copy}')  # pylint:disable=C0103
                best_alt_prob = np.max(ap, axis=1)
                ref_prob = 1 - np.sum(ap, axis=1)
                best_copy_prob[copy] = np.max(
                    (best_alt_prob, ref_prob),
                    axis=0
                )
            probs[variant_count, :] = 10000 * \
                np.multiply(best_copy_prob[1], best_copy_prob[2])


def _zarr_store():
    return zarr.LMDBStore(zarr_store_loc)


# zpath stands for zarr path
def _phased_hardcall_probs_zpath(chrom):
    return f'chr{chrom}/phased_hardcall_probs'


def _phased_hardcalls_zpath(chrom):
    return f'chr{chrom}/phased_hardcalls'


# TODO how am I going to parallelize this?
# a different database per chromosome for additional write parallelization?
# or just one database across all?
def main():
    print("Creating store")
    sys.stdout.flush()
    create_phased_hardcall_store(487409, 1000, 1)
    ukb = os.environ['UKB']
    print("Staring to read vcf")
    sys.stdout.flush()
    read_phased_hardcalls_from_vcf(
        f'{ukb}/str_imputed/runs/first_pass/vcfs/strs_only/chr1.vcf.gz',
        1
    )


if __name__ == '__main__':
    main()
