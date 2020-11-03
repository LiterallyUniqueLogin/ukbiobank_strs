# pylint: disable=C0103,C0114,C0116
import os
import os.path
import subprocess as sp

# Use Dask to smooth the PBS submission process
import dask.bag
import dask.distributed
import dask_jobqueue

ukb = os.environ['UKB']

def download_item(item_info):  # noqa: D103
    ukb, chrom, bgen_dir = item_info
    command = f"""
    cd {bgen_dir} ;
    {ukb}/ukb_utilities/ukbgene \
            imp \
            -a../main_dataset/raw_data/k29170.key \
            -c{chrom}
    """
    try :
        output = sp.run(command,
                        shell=True,
                        check=True,
                        stdout=sp.PIPE,
                        stderr=sp.PIPE)
        print(output.stdout.decode(), flush=True)
    except sp.CalledProcessError as cpe:
        print("----stdout---")
        print(cpe.stdout.decode())
        print("----stderr---")
        print(cpe.stderr.decode(), flush=True)
        raise cpe

def main():  # noqa: D103

    bgen_dir = f'{ukb}/array_imputed'
    output_dir = f'{bgen_dir}/output'

    assert os.path.exists(bgen_dir)
    assert os.path.exists(output_dir)

    # why not retry? nothing else I can do
    dask.config.set({'distributed.scheduler.allowed-failures': 99})

    # Maximum of 10 concurrent downloads per application
    # See here: https://biobank.ctsu.ox.ac.uk/showcase/refer.cgi?id=644
    client = dask.distributed.Client(
        n_workers=10,
        local_directory="/oasis/tscc/scratch/jmargoli"
    )

    jobs = []
    # calculate number of download batches
    for chrom in range(1,23):
        jobs.append((ukb, chrom, bgen_dir))

    print(f"Number of jobs queued: {len(jobs)}", flush = True)

    bag = dask.bag.from_sequence(jobs)
    downloads = bag.map(download_item)

    client.compute(downloads, retries=99).result() # wait for the result so
    # so that all tasks complete

if __name__ == "__main__":
    main()
