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
    ukb, sample_ID, field_ID, vcf_dir = item_info
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
        print(output.stdout.decode(), flush=True)
    except sp.CalledProcessError as cpe:
        print("----stdout---")
        print(cpe.stdout.decode())
        print("----stderr---")
        print(cpe.stderr.decode(), flush=True)
        raise cpe

def main():  # noqa: D103

    tscc_vcf_dir = f'{ukb}/../../resources/datasets/ukbiobank/exome/fe_crams'
    vcf_dir = tscc_vcf_dir
    bulk_floc = f'{ukb}/exome/fe_cram.bulk'

    assert os.path.exists(vcf_dir)
    assert os.path.exists(bulk_floc)

    current_files = set(os.listdir(tscc_vcf_dir))

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
    with open(bulk_floc) as bulk_file:
        for line in bulk_file:
            sample_ID, field_ID = line.split()
            if field_ID == '23163_0_0':
                suffix = 'cram'
            elif field_ID == '23164_0_0':
                suffix = 'cram.crai'
            file_name = f"{sample_ID}_{field_ID}.{suffix}"
            if file_name in current_files:
                continue
            jobs.append((ukb, sample_ID, field_ID, vcf_dir))

    print(f"Number of jobs queued: {len(jobs)}", flush = True)

    bag = dask.bag.from_sequence(jobs)
    downloads = bag.map(download_item)

    client.compute(downloads, retries=99).result() # wait for the result so
    # so that all tasks complete

if __name__ == "__main__":
    main()
