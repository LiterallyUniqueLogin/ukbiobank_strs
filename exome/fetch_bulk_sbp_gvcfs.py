# pylint: disable=C0103,C0114,C0116
import math
import os
import subprocess as sp

# Use Dask to smooth the PBS submission process
import dask
import dask.distributed
import dask_jobqueue

BATCH_SIZE = 1000


@dask.delayed
def download_items(batch):  # noqa: D103
    command = f"""
    cd {vcf_dir}
    ../../ukb_utilities/ukbfetch \
            -b../spb_gvcfs.bulk \
            -a../../main_dataset/k41414.key \
            -v \
            -of{batch} \
            -s{batch * BATCH_SIZE} \
            -m{BATCH_SIZE}
    """

    try:
        output = sp.run(command,
                        shell=True,
                        check=True,
                        stdout=sp.PIPE,
                        stderr=sp.PIPE)
        print(output.stdout.decode())
        return None
    except sp.CalledProcessError as e:
        return (e.returncode, e.cmd, e.stdout, e.stderr)


def main():  # noqa: D103
    output_dir = f'{vcf_dir}/output'
    cluster = dask_jobqueue.PBSCluster(
        name="STR_subset",
        walltime="24:00:00",
        log_directory=output_dir
    )
    cluster.scale(20)
    client = dask.distributed.Client(cluster)

    # calculate number of download batches
    with open(f'{vcf_dir}/../spb_gvcfs.bulk') as bulk_file:
        nitems = len(bulk_file.readlines())

    jobs = []
    for batch in range(0, math.ceil(nitems / BATCH_SIZE)):
        jobs.append(download_items(batch))

    futures = client.compute(jobs)

    with open(f"{output_dir}/results.txt", 'w') as results_file:
        for batch, future in enumerate(futures):
            result = future.result()
            if result is None:
                results_file.write(f"batch {batch} succeeded\n")
            else:
                results_file.write(
                    f"batch {batch} failed. Error: {result}\n\n"
                )


if __name__ == "__main__":
    ukb = os.environ['UKB']
    vcf_dir = f'{ukb}/exome/vcfs'
    main()
