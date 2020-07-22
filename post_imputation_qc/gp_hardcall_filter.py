import argparse
import os
import os.path
import subprocess as sp
import sys
import time

import cyvcf2
import dask.distributed
import dask_jobqueue
import numpy as np
import numpy.lib.arraysetops as setops
import zarr

# Use a LMDB store
# This seems to suggest that LMDB is faster than simply using the filesystem
# https://groups.google.com/forum/#!topic/caffe-users/WqBfwUJ98c0
# Confirm this

# Zarr's default compression algorithm, blosc, uses 8 threads internally. Make
# sure this is respected when submitting this program to the cluster

# create array in Zarr big enough to store all the data we want

ukb = os.environ['UKB']

zarr_store_name = 'filtering.lmdb.zarr'
zarr_store_size = 2**38  # 274GB
max_n_alleles, n_allele_dtype = 2**8 - 1, 'uint8'


def get_sample_bitmap(vcf_floc, sample_floc):
    """
    Return
    ------
    A 1-D array bit-array, one entry per sample in the input vcf
    (same ordering)
    with each bit indicating the presence or absence of that sample
    in the samples file (1 indicates present)
    """
    with open(sample_floc) as sample_file:
        samples_iter = iter(sample_file)
        first_line = next(samples_iter)
        assert first_line.startswith('ID')
        included_samples = []
        for line in samples_iter:
            sample = line.split()[0]
            included_samples.append(f'{sample}_{sample}')

    included_samples = np.array(included_samples)
    vcf = cyvcf2.VCF(vcf_floc)

    return setops.isin(vcf.samples, included_samples)


def count_variants_for_chrom(chrom, vcf_floc):
    # use bcftools index -n to count variants
    command = f"""
source ~/.bashrc
conda activate bcftools
bcftools index -n {vcf_floc}
conda deactivate
"""
    output = sp.run(command,
                    shell=True,
                    check=True,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE)
    return int(output.stdout.decode())


class mem_buf_writing_zarr_array:
    def __init__(self, zarr_array, nrows_per_chunk, ncols, dtype):
        self.zarr_array = zarr_array
        self.nrows_per_chunk = nrows_per_chunk
        self.ncols = ncols
        self.dtype = dtype
        self.row_chunk_idx = 0
        self.row_idx = 0
        self._reset_mem_buf()
        
    def _reset_mem_buf(self):
        self.mem_buf = np.full(
            (self.nerows_per_chunk, self.ncols),
            np.nan,
            dtype=self.dtype
        )

    def fill_row(self, row):
        self.mem_buf[self.row_idx, :] = row
        self.row_idx += 1
        if self.row_idx == self.nrow_per_chunk:
            self.row_idx = 0
            self.zarr_array[(self.row_chunk_idx * self.nrows_per_chunk): \
                            ((self.row_chunk_idx + 1) * self.nrows_per_chunk),
                            :] = self.mem_buf
            self.row_chunk_idx += 1
            self._reset_mem_buf()


def create_phased_hardcall_store(filtering_run_dir,
                                 n_variants,
                                 n_samples,
                                 chrom):
    with _zarr_store(filtering_run_dir, chrom) as store:
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


def read_phased_hardcalls_from_vcf(filtering_run_dir,
                                   vcf_floc,
                                   chrom,
                                   sample_bitmap):
    with _zarr_store(filtering_run_dir, chrom) as store:
        gts_array = zarr.open_array(
            store=store,
            path=_phased_hardcalls_zpath(chrom)
        )
        gts = mem_buf_writing_zarr_array(gts_array)
        probs_array = zarr.open_array(
            store=store,
            path=_phased_hardcall_probs_zpath(chrom)
        )
        probs = mem_buf_writing_zarr_array(probs_array)

        vcf = cyvcf2.VCF(vcf_floc)

        start = time.time()
        prev = time.time()
        variant_count = 0
        allele_count = 0
        prev_allele_count = 0
        logging_interval = 20  # in variants
        for variant in vcf:
            if len(variant.ALT) > max_n_alleles - 1:
                raise OverflowError((
                    f"Variant {variant.POS} had {len(variant.ALT)}alternate "
                    f"alleles but only had space for {max_n_alleles - 1} "
                    f"alternate alleles"
                ))

            # set the genotypes
            gts.fill_row(variant.genotype.array()[sample_bitmap, :2])

            # set the probabilities
            best_copy_prob = {}
            for copy in [1, 2]:
                ap = variant.format(f'AP{copy}')[sample_bitmap, :]
                best_alt_prob = np.max(ap, axis=1)
                ref_prob = 1 - np.sum(ap, axis=1)
                best_copy_prob[copy] = np.max(
                    (best_alt_prob, ref_prob),
                    axis=0
                )
            probs.fill_row(10000 * np.multiply(best_copy_prob[1],
                                               best_copy_prob[2]))

            variant_count += 1
            prev_allele_count += len(variant.ALT) + 1
            allele_count += len(variant.ALT) + 1

            if variant_count % logging_interval == 0:
                now = time.time()
                interval_time = now - prev
                interval_time_per_variant = interval_time / logging_interval
                interval_time_per_allele = interval_time / prev_allele_count
                total_time = now - start
                total_time_per_variant = total_time / variant_count
                total_time_per_allele = total_time / allele_count
                print(f"Finished variant {variant_count} at POS {variant.POS}")
                print(f"Time per variant last {logging_interval} variants: "
                      f"{interval_time_per_variant:.2f}sec")
                print(f"Time per allele last {logging_interval} "
                      f"variants: {interval_time_per_allele:.2f}sec")
                print(f"Time per variant over all variants: "
                      f"{total_time_per_variant:.2f}sec")
                print(f"Time per allele over all variants: "
                      f"{total_time_per_allele:.2f}sec")
                sys.stdout.flush()
                prev_allele_count = 0
                prev = time.time()


def _zarr_store(filtering_run_dir, chrom):
    os.makedirs(f"{filtering_run_dir}/{zarr_store_name}",
                exist_ok=True)
    return zarr.LMDBStore(
        f"{filtering_run_dir}/{zarr_store_name}/chr{chrom}",
        map_size=zarr_store_size
        # writemap = True?
    )


# zpath stands for zarr path
def _phased_hardcall_probs_zpath(chrom):
    return f'chr{chrom}/phased_hardcall_probs'


def _phased_hardcalls_zpath(chrom):
    return f'chr{chrom}/phased_hardcalls'


@dask.delayed
def load_phased_hardcalls_into_zarr(chrom,
                                    filtering_run_dir,
                                    vcf_floc,
                                    samples_floc):
    print(f"Processing chromosome {chrom}")
    print("Collecting num variants, samples")
    sys.stdout.flush()
    sample_bitmap = get_sample_bitmap(vcf_floc, samples_floc)
    n_samples = np.sum(sample_bitmap)
    n_variants = count_variants_for_chrom(chrom, vcf_floc)
    print(f"n_variants {n_variants}, n_samples {n_samples}")
    print("Creating store")
    sys.stdout.flush()
    create_phased_hardcall_store(filtering_run_dir,
                                 n_variants,
                                 n_samples,
                                 chrom)
    print("Staring to read vcf")
    sys.stdout.flush()
    read_phased_hardcalls_from_vcf(filtering_run_dir,
                                   vcf_floc,
                                   chrom,
                                   sample_bitmap)
    print("Done")


def main():  # noqa: D103
    parser = argparse.ArgumentParser()
    parser.add_argument('filtering_run_name')
    parser.add_argument('imputation_run_name')
    parser.add_argument('samples_file')
    args = parser.parse_args()

    filtering_run_name = args.filtering_run_name
    imputation_run_name = args.imputation_run_name
    samples_file = args.samples_file

    filtering_run_dir = f'{ukb}/post_imputation_qc/runs/{filtering_run_name}'
    output_dir = f'{filtering_run_dir}/output'
    os.makedirs(output_dir, exist_ok=True)

    with open(f'{filtering_run_dir}/README', 'w') as readme:
        readme.write(f'imputation_run_name={imputation_run_name}\n'
                     f'samples_file={samples_file}\n')

    cluster = dask_jobqueue.PBSCluster(
        name=f"filtering_{filtering_run_name}",
        walltime="48:00:00",
        log_directory=output_dir
    )
    cluster.scale(22)
    client = dask.distributed.Client(cluster)

    jobs = []
    for chrom in range(1, 23):
        vcf_floc = (f'{ukb}/str_imputed/runs/{imputation_run_name}/'
                    f'vcfs/strs_only/chr{chrom}.vcf.gz')
        jobs.append(load_phased_hardcalls_into_zarr(
            chrom,
            filtering_run_dir,
            vcf_floc,
            samples_file
        ))

    futures = client.compute(jobs)

    with open(f"{output_dir}/results.txt", 'w') as results_file:
        for chrom, future in enumerate(futures):
            chrom += 1
            err = future.exception()
            if err is None:
                results_file.write(f"chrom {chrom} succeeded\n\n")
            else:
                results_file.write(
                    f"chrom {chrom} failed. Error: {err}\n\n"
                )


if __name__ == '__main__':
    main()
