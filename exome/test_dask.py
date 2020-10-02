# Use Dask to smooth the PBS submission process
import dask
import dask.distributed
import dask_jobqueue
import dask.delayed

import os

@dask.delayed
def fail():
   raise ValueError()

def main():  # noqa: D103
    output_dir = f'{os.environ["UKB"]}/exome/test_output'

    cluster = dask_jobqueue.PBSCluster(
        name="test",
        walltime="4:00:00",
        log_directory=output_dir,
        queue="condo"
    )
    # Maximum of 10 concurrent downloads per application
    # See here: https://biobank.ctsu.ox.ac.uk/showcase/refer.cgi?id=644
    cluster.adapt(minimum_jobs=10, maximum_jobs=10)
    client = dask.distributed.Client(cluster)

    jobs = set()
    jobs.add(fail())

    futures = client.compute(jobs, retries=1)

    for future in futures:
        future.result() # block till code is done executing
        results_file.write(f"{future.key} succeeded\n")

if __name__ == "__main__":
    main()
