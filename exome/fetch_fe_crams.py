# pylint: disable=C0103,C0114,C0116
import datetime
import os
import os.path
import subprocess as sp

# Use Dask to smooth the PBS submission process
import dask
import dask.distributed
import dask_jobqueue
import dask.delayed

ukb = os.environ['UKB']

@dask.delayed
def download_item(ukb, sample_ID, field_ID, vcf_dir):  # noqa: D103
    command = f"""
    cd {vcf_dir} ;
    {ukb}/ukb_utilities/ukbfetch \
            -ak41414.key \
            -v \
            -of{sample_ID}_{field_ID} \
            -e{sample_ID} \
            -d{field_ID}
    """
    try :
        output = sp.run(command,
                        shell=True,
                        check=True,
                        stdout=sp.PIPE,
                        stderr=sp.PIPE)
        print(output.stdout.decode())
    except sp.CalledProcessError as cpe:
        print("----stdout---")
        print(cpe.stdout)
        print("----stderr---")
        print(cpe.stderr)

def main():  # noqa: D103
    vcf_dir = f'{ukb}/exome/fe_crams'
    output_dir = f'{vcf_dir}_output'
    bulk_floc = f'{ukb}/exome/fe_cram.bulk'

    assert os.path.exists(vcf_dir)
    assert os.path.exists(bulk_floc)

    current_files = set(os.listdir(vcf_dir))
    cluster = dask_jobqueue.PBSCluster(
        name="UKB_fe_cram_download",
        walltime="4:00:00",
        log_directory=output_dir,
        queue="condo"
    )
    # Maximum of 10 concurrent downloads per application
    # See here: https://biobank.ctsu.ox.ac.uk/showcase/refer.cgi?id=644
    cluster.adapt(minimum_jobs=10, maximum_jobs=10)
    client = dask.distributed.Client(cluster)

    jobs = set()
    # calculate number of download batches
    with open(bulk_floc) as bulk_file:
        for line in bulk_file:
            sample_ID, field_ID = line.split()
            if field_ID == '23163_0_0':
                suffix = 'cram'
            elif field_ID == '23164_0_0':
                suffix = 'crai'
            file_name = f"{sample_ID}_{field_ID}.{suffix}"
            if file_name in current_files:
                continue
            jobs.add(download_item(
                ukb, sample_ID, field_ID, vcf_dir,
                dask_key_name=f'download_item-{sample_ID}-{field_ID}'
            ))

    print(f"Number of jobs queued: {len(jobs)}")

    futures = client.compute(jobs, retries=1)

    now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    with open(f"{output_dir}/results_{now}.txt", 'w') as results_file:
        for future in futures:
            future.result() # block till code is done executing
            results_file.write(f"{future.key} succeeded\n")

if __name__ == "__main__":
    main()
