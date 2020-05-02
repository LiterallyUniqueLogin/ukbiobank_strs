import os
import os.path

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
    if os.path.exists(zarr_store_loc):
        raise FileExistsError(zarr_store_loc)

    # hardcalls are calls which take the genotype with the highest probability
    # and ignore the probabilities of the rest of the genotypes.
    # phased hardcalls mean that '0|1' and '1|0' are different genotypes and
    # so their probabilities are not combined. These are the GTs that Beagle
    # outputs. Note that it is probably more realistic to work with unphased
    # hardcalls (i.e. p('0/1') which equals p('0|1') + p('1|0')) but
    # calculating those would be more compute intensive and the resulting
    # difference is unlikely to matter

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
        store=_zarr_store(),
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
        store=_zarr_store(),
        path=_phased_hardcalls_zpath(chrom)
    )


def read_phased_hardcalls_from_vcf(vcf_floc, chrom):
    gts = zarr.open_array(
        store=_zarr_store(),
        path=_phased_hardcalls_zpath(chrom)
    )
    probs = zarr.open_array(
        store=_zarr_store(),
        path=_phased_hardcall_probs_zpath(chrom)
    )
    vcf = cyvcf2.VCF(vcf_floc)

    count = -1
    for variant in vcf:
        count += 1

        if len(variant.ALT) > max_n_alleles - 1:
            raise OverflowError((
                f"Variant {variant.POS} had {len(variant.ALT)}alternate "
                f"alleles but only had space for {max_n_alleles - 1} "
                f"alternate alleles"
            ))

        # set the genotypes
        gts[count, :, :] = variant.genotype.array()[:, :2]

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
        probs[count, :] = 10000 * \
            np.multiply(best_copy_prob[1], best_copy_prob[2])
        break


_ZARR_STORE = None


def _zarr_store():
    global _ZARR_STORE  # pylint:disable=W0603
    if _ZARR_STORE is None:
        _ZARR_STORE = zarr.LMDBStore(zarr_store_loc)
    return _ZARR_STORE


def _close_store():
    _ZARR_STORE.close()


# zpath stands for zarr path
def _phased_hardcall_probs_zpath(chrom):
    return f'chr{chrom}/phased_hardcall_probs'


def _phased_hardcalls_zpath(chrom):
    return f'chr{chrom}/phased_hardcalls'


if __name__ == '__main__':
    create_phased_hardcall_store(487409, 1000, 1)
    ukb = os.environ['UKB']
    read_phased_hardcalls_from_vcf(
        f'{ukb}/str_imputed/runs/first_pass/vcfs/strs_only/chr1.vcf.gz',
        1
    )
