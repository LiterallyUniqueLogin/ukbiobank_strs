# pylint: disable=C0103,C0114,C0116
import os
import subprocess as sp

# Use Dask to smooth the PBS submission process
import dask
import dask.distributed
import dask_jobqueue


@dask.delayed
def download_item(sample_ID, field_ID):  # noqa: D103
    command = f"""
    cd {vcf_dir}
    ../../ukb_utilities/ukbfetch \
            -a../../main_dataset/k41414.key \
            -v \
            -of{sample_ID}_{field_ID} \
            -e{sample_ID} \
            -d{field_ID}
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
    current_files = set(os.listdir(vcf_dir))
    cluster = dask_jobqueue.PBSCluster(
        name="UKB_gVCF_download",
        walltime="4:00:00",
        log_directory=output_dir
    )
    # Maximum of 10 concurrent downloads per application
    # See here: https://biobank.ctsu.ox.ac.uk/showcase/refer.cgi?id=644
    cluster.adapt(minimum_jobs=10, maximum_jobs=10)
    client = dask.distributed.Client(cluster)

    jobs = set()
    # calculate number of download batches
    with open(f'{vcf_dir}/../spb_gvcfs.bulk') as bulk_file:
        for line in bulk_file:
            sample_ID, field_ID = line.split()
            if field_ID == '23176_0_0':
                suffix = 'gz'
            elif field_ID == '23177_0_0':
                suffix = 'tbi'
            file_name = f"{sample_ID}_{field_ID}.{suffix}"
            if file_name in current_files:
                continue
            jobs.add(download_item(sample_ID, field_ID))

    print(f"Number of jobs queued: {len(jobs)}")

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
    output_dir = f'{vcf_dir}/output'
    main()
