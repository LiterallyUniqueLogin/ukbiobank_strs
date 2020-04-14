# pylint: disable=C0103,C0114,C0116
import argparse
import os
import subprocess as sp

# Use Dask to smooth the PBS submission process
import dask
import dask.distributed
import dask_jobqueue


@dask.delayed
def subset_strs(chrom):  # noqa: D103
    subset_vcf = f'{subset_dir}/chr{chrom}.vcf.gz'
    ref_loc = (f'{dataset_loc}/1000Genomes/hipstr_calls_sample_subset/'
               f'eur/chr{chrom}/hipstr.chr{chrom}.eur.vcf.gz')
    command = f"""
    source ~/.bashrc
    conda activate bcftools
    bcftools view \
        -O z \
        -o {subset_dir}/chr{chrom}.vcf.gz \
        -r {ukb}/snpstr/str_loci.txt \
        -i ID=@{ukb}/snpstr/str_ids.txt \
        {ukb}/str_imputed/runs/{args.run_name}/vcfs/chr{chrom}.vcf.gz
    tabix {subset_vcf}
    last_subset_id=$(bcftools query -f '%ID\n' {subset_vcf} | tail -n 1)
    echo "last_subset_id $last_subset_id"
    all_last_ids=$(bcftools query -f '%ID\n' {ref_loc} | tail -n 300)
    if ! echo "$all_last_ids" | grep $last_subset_id ; then
        echo "Last subset id $last_subset_id is not towards the end of \
the list of STR ids" >&2
        exit 1
    fi
    conda deactivate
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
    os.makedirs(output_dir, exist_ok=True)

    cluster = dask_jobqueue.PBSCluster(
        name="STR_subset",
        walltime="48:00:00",
        log_directory=output_dir
    )
    cluster.scale(22)
    client = dask.distributed.Client(cluster)

    jobs = []
    for chrom in range(1, 23):
        jobs.append(subset_strs(chrom))

    futures = client.compute(jobs)

    with open(f"{output_dir}/results.txt", 'w') as results_file:
        for chrom, future in enumerate(futures):
            chrom += 1
            result = future.result()
            if result is None:
                results_file.write(f"chrom {chrom} succeeded\n\n")
            else:
                results_file.write(
                    f"chrom {chrom} failed. Error: {result}\n\n"
                )


if __name__ == "__main__":
    ukb = os.environ['UKB']
    dataset_loc = os.environ['DATASETS']

    parser = argparse.ArgumentParser()
    parser.add_argument("run_name")
    args = parser.parse_args()
    subset_dir = f"{ukb}/str_imputed/runs/{args.run_name}/vcfs/strs_only"
    output_dir = f"{subset_dir}/output"
    main()
