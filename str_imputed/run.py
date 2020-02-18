import argparse
import datetime
import os
import os.path
import subprocess as sp
import sys

sys.path.insert(0, os.environ['UKB'] + "/str_imputed/scripts")

import check_beagle_output_samples
import check_final_output_samples
import check_merge_output_samples
import launch_impute
import launch_merge_within_region

parser = argparse.ArgumentParser(
    description=("Continues the imputation process by assessing which jobs "
                 "need to be run and launching them. Information about the "
                 "current status of the imputation, issues that need to be "
                 "resolved, and next steps, to "
                 "$UKB/str_imputed/runs/<run_name>/status_<date>.txt"))
parser.add_argument(
    "run_name",
    help=("The name of this run. Output files will be put in "
          "$UKB/str_imputed/runs/<run_name>"))
parser.add_argument(
    "--pfile_directory",
    help=("The directory containing the pfiles of the dataset to be imputed "
          "(the pfiles must be named chr1 ... chr22). Only necessary the "
          "first time running this command"))
parser.add_argument(
    "--sample_file",
    help=("the sample file with the list of samples, see "
          "$UKB/microarray/*.sample for an example. Only necessary the first "
          "time running this command"))
parser.add_argument(
    "--readme",
    help=("the description of this run (e.g. what filters were use to create "
          "the input .sample file and pfiles, etc.) Only necessary the first "
          "time running this command"))
parser.add_argument(
    "--debug",
    action='store_true',
    help="Skip collecting error files because that is slow",
    default=False)

args = parser.parse_args()
run_name = args.run_name
pfile_directory = args.pfile_directory
sample_file = args.sample_file
readme = args.readme

def error(msg):
    print(msg, file=sys.stderr)
    exit(-1)

if 'TMPDIR' not in os.environ:
    error("Must set the TMPDIR environment variable")

run_dir = os.environ['UKB'] + "/str_imputed/runs/" + run_name
os.makedirs(run_dir, exist_ok=True)
os.makedirs(f"{run_dir}/batches", exist_ok=True)
os.makedirs(f"{run_dir}/batches/output", exist_ok=True)
os.makedirs(f"{run_dir}/batches/old", exist_ok=True)
os.makedirs(f"{run_dir}/vcfs", exist_ok=True)

# Load num samples
info_file = f'{run_dir}/info.txt'
if not os.path.exists(info_file):
    if not pfile_directory or not sample_file:
        error("Options --pfile_directory and --sample_directory are mandatory "
              "on first run")

    tmpdir = os.environ['TMPDIR']
    os.makedirs(f'{run_dir}', exist_ok=True)
    sp.run(f"plink2 --pfile {pfile_directory}/chr1 --keep {sample_file} "
           f"--write-samples --out {tmpdir}/{run_name}",
           stdout=sp.PIPE, shell=True, check=True)
    stdout = sp.run(f"wc -l {tmpdir}/{run_name}.id", check=True,
                    shell=True, stdout=sp.PIPE, stderr=sp.PIPE).stdout
    num_samples = int(stdout.decode().split()[0]) - 1
    if num_samples < 50000:
        error(f"Error: found too few samples ({num_samples}), something went "
              "wrong during pfile loading")

    with open(info_file, 'w') as info:
        info.write(f"pfile_directory:{pfile_directory}\n"
                   f"sample_file:{sample_file}\n"
                   f"num_samples:{num_samples}\n")
else:
    with open(info_file) as info:
        fields = {}
        for line in info:
            key, value = line.split(':')
            fields[key.strip()] = value.strip()
        pfile_directory = fields['pfile_directory']
        sample_file = fields['sample_file']
        num_samples = int(fields['num_samples'])

scripts_dir = os.environ['UKB'] + "/str_imputed/scripts"
final_files = os.listdir(run_dir + '/vcfs')
batch_files = os.listdir(run_dir + '/batches')
output_files = os.listdir(run_dir + '/batches/output')

if not args.debug:
    print("Collecting errors")
    errors_output = sp.run(f'{scripts_dir}/check_for_pbs_errors.sh {run_name} '
                           ' | grep "Error in file"',
                           check=True,
                           shell=True,
                           universal_newlines=True,
                           stdout=sp.PIPE)
    error_files = set()
    for line in errors_output.stdout.split('\n'):
        if line.strip() == "":
            continue
        error_files.add(line.split(' ')[3])
else:
    error_files = set()

job_names = set()
job_ids_output = sp.run('qselect -u $(whoami)',
                        shell=True, universal_newlines=True,
                        stdout=sp.PIPE, check=True).stdout.split('\n')
job_ids = " ".join(job_id.split('.')[0] for job_id in job_ids_output)
qstat_output = sp.run(f'qstat -f {job_ids} | grep Job_Name ',
                      shell=True, universal_newlines=True,
                      stdout=sp.PIPE, check=True).stdout
for line in qstat_output.split("\n"):
    if line.strip() == "":
        continue
    job_names.add(line.split("=")[1].strip())

def check_impute(status,
                 batch_files,
                 job_names,
                 error_files,
                 output_files,
                 chrs_to_check):

    status.write("---Step 1: Imputation---\n")
    keep_checking = set()
    done_with_impute_file = f"{run_dir}/done_with_1_impute.txt"
    if os.path.exists(done_with_impute_file):
        with open(done_with_impute_file) as done_with_impute:
            keep_checking = \
                keep_checking.union(int(chrom) for chrom in
                                    done_with_impute.read().split('\n')
                                    if chrom != "")

    for chrom in chrs_to_check:
        if chrom in keep_checking:
            continue

        print(f"Checking chrom {chrom} impute progress")

        # run launch_impute which handles errors itself
        currently_running_jobs = set()
        for sample_idx in range(1, num_samples + 1, 1000):
            job_name = f"impute_STRs_{chrom}_{sample_idx}"
            if job_name in job_names:
                currently_running_jobs.add(job_name)
        ret_code = launch_impute.do_impute(
            run_name, pfile_directory, sample_file,
            chrom, readme, ignore_jobs=currently_running_jobs)
        if type(ret_code) == str:
            status.write(f"Launching imputation jobs for chromosome {chrom} "
                         f"resulted in the error {ret_code}\n")
            continue

        if ret_code is False:
            status.write(f"Chromosome {chrom} is running imputations with "
                         f"names beginning with impute_STRs_{chrom}\n")
            continue

        ret_code = check_beagle_output_samples.do_check(
            run_name, sample_file, chrom, command_line=False)
        if type(ret_code) == str:
            status.write(f"Chromosome {chrom} was imputed with the wrong "
                         f"samples! Please fix this manually. Specific error: "
                         f"{ret_code}\n")
            continue

        # check imputation variants are correct
        job_name = f"check_beagle_output_variants_{chrom}"
        if job_name in job_names:
            status.write(f"Chrom {chrom} is running the beagle output check "
                         f"variants job, with name {job_name}\n")
            continue

        error_file = name_from_tag(job_name, error_files)
        if error_file:
            status.write(f"Chrom {chrom} failed the beagle output variants "
                         f"check. Errors detailed in file f{error_file}\n")
            continue

        output_file = name_from_tag(f"{job_name}.o", output_files)
        if output_file is None:
            status.write(f"Chrom {chrom} is starting up the beagle output "
                         f"check variants job, with name {job_name}\n")
            sp.run(f"{scripts_dir}/launch_check_beagle_output_variants.sh "
                   f"{run_name} {chrom}",
                   shell=True, stdout=sp.PIPE, check=True)
        else:
            output_file = f"{run_dir}/batches/output/{output_file}"
            with open(output_file) as output:
                if "Success" not in output.read():
                    status.write(f"Chrom {chrom} was tested for imputation "
                                 f"output variants but the testing didn't "
                                 f"error or return success. Please diagnose "
                                 f"manually. See file {output_file}\n")
                    continue

            keep_checking.add(chrom)

    with open(done_with_impute_file, 'w') as done_with_impute:
        first = True
        for chrom in keep_checking:
            if not first:
                done_with_impute.write("\n")
            first = False
            done_with_impute.write(str(chrom))

    return keep_checking

def check_regional_merge(status,
                         batch_files,
                         job_names,
                         error_files,
                         output_files,
                         chrs_to_check):
    status.write("---Step 2: Regional Merge---\n")
    keep_checking = set()
    done_with_regional_merge_file = f"{run_dir}/done_with_2_regional_merge.txt"
    if os.path.exists(done_with_regional_merge_file):
        with open(done_with_regional_merge_file) as done_with_regional_merge:
            keep_checking = \
                keep_checking.union(
                    int(chrom) for chrom in
                    done_with_regional_merge.read().split('\n')
                    if chrom != "")

    for chrom in chrs_to_check:
        if chrom in keep_checking:
            continue

        print(f"Checking chrom {chrom} regional merge progress")

        keep_checking_this_chrom = True

        sample_failures = check_merge_output_samples.do_check(
            run_name, chrom, command_line=False)

        currently_running_jobs = False
        errored_jobs = False
        final_sample_failures = set()
        new_jobs = set()
        check_num_variants = False
        check_num_variants_errors = False
        new_checks = set()

        for pos in range(1, 250000000, 5000000):
            job_name = f"merge_within_region_{chrom}_{pos}"
            if job_name in job_names:
                currently_running_jobs = True
                keep_checking_this_chrom = False
                continue

            temp_error_file = name_from_tag(job_name, error_files)
            if temp_error_file is not None:
                error_file = temp_error_file
                errored_jobs = True
                keep_checking_this_chrom = False
                continue

            batch_file = name_from_tag(f"chr{chrom}_pos_{pos}_", batch_files)
            if batch_file is None:
                new_jobs.add(str(pos))
                keep_checking_this_chrom = False
                continue

            if pos in sample_failures:
                final_sample_failures.add(pos)
                keep_checking_this_chrom = False
                continue

            job_name = f"check_merge_output_num_variants_{chrom}_{pos}"
            if job_name in job_names:
                keep_checking_this_chrom = False
                check_num_variants = True
                continue

            temp_check_error_file = name_from_tag(job_name, error_files)
            if temp_check_error_file is not None:
                check_error_file = temp_check_error_file
                keep_checking_this_chrom = False
                check_num_variants_errors = True
                continue

            output_file = name_from_tag(f"{job_name}.o", output_files)
            if output_file is None:
                new_checks.add(str(pos))
                keep_checking_this_chrom = False
            else:
                with open(f"{run_dir}/batches/output/{output_file}") as output:
                    if 'Success' not in output.read():
                        status.write(f"Chromosome {chrom} check regional "
                                     f"merge num varaints jobs finished "
                                     f"without either error or success. "
                                     f"Please debug manually. See file "
                                     f"{output_file}\n")
                        keep_checking_this_chrom = False
                # this position is fully successful, don't need to do anything

        if currently_running_jobs:
            status.write(f"Chromosome {chrom} is running regional merges\n")
        if errored_jobs:
            status.write(f"Chromosome {chrom} failed the regional merge in "
                         f"at least one region. Errors detailed in files "
                         f"named similarly to f{error_file}\n")

        if len(new_jobs) > 0:
            msg = launch_merge_within_region.do_launch(
                run_name, chrom, new_jobs, command_line=False)
            if msg is not None:
                status.write(f"Chromosome {chrom} failed to launch regional "
                             f"merge jobs. Please debug manually. See error "
                             f"{msg}\n")
            else:
                status.write(f"Chromsome {chrom} is launching new regional "
                             f"merge jobs\n")

        if len(final_sample_failures) != 0:
            status.write(f"Chromosome {chrom} was merged with the wrong "
                         f"samples in regions {final_sample_failures}! Please "
                         f"fix this manually.\n")

        if check_num_variants:
            status.write(f"Chrom {chrom} is running merge output check "
                         f"variants job(s), with names beginning with "
                         f"check_merge_output_num_varaints_{chrom}\n")
        if check_num_variants_errors:
            status.write(f"Chromosome {chrom} failed the regional merge num "
                         f"variants check in at least one region. Errors "
                         f"detailed in files named similarly to "
                         f"{check_error_file}\n")

        if len(new_checks) > 0:
            regions = " ".join(list(new_checks))
            sp.run(f"{scripts_dir}/launch_check_merge_output_num_variants.sh "
                   f"{run_name} {chrom} {regions}",
                   shell=True, stdout=sp.PIPE, check=True)
            status.write(f"Chromosome {chrom} is running new regional merge "
                         f"num variants checks. Job names similar to "
                         f"check_merge_output_num_variants_{chrom}\n")

        if keep_checking_this_chrom:
            keep_checking.add(chrom)

    with open(done_with_regional_merge_file, 'w') as done_with_regional_merge:
        first = True
        for chrom in keep_checking:
            if not first:
                done_with_regional_merge.write("\n")
            first = False
            done_with_regional_merge.write(str(chrom))

    return keep_checking

def check_concat(status,
                 final_files,
                 job_names,
                 error_files,
                 output_files,
                 chrs_to_check):
    status.write("---Step 3: Concat---\n")
    for chrom in chrs_to_check:

        print(f"Checking chrom {chrom} concat progress")

        job_name = f"concat_chr_regions_{chrom}"
        if job_name in job_names:
            status.write(f"The concat job for chrom {chrom} is currently "
                         f"running with name {job_name}\n")
            continue

        error_file = name_from_tag(job_name, error_files)
        if error_file:
            status.write(f"Concat for chrom {chrom} failed. Errors detailed "
                         f"in file f{error_file}\n")
            continue

        tbi_file = f'chr{chrom}.vcf.gz.tbi'
        if f'chr{chrom}.vcf.gz' not in final_files or \
                tbi_file not in final_files:
            status.write(f"Launching the concat job for chrom {chrom}, job "
                         f"name concat_chr_regions_{chrom}\n")
            sp.run(f"{scripts_dir}/launch_concat_chr_regions.sh "
                   f"{run_name} {chrom}",
                   shell=True, stdout=sp.PIPE, check=True)
            continue

        if os.stat(run_dir + "/vcfs/" + tbi_file).st_size < 2e4:
            status.write(f"Relaunching the concat job for chrom {chrom} - "
                         f"the tabix index file produced is smaller than "
                         f"anticipated, inticating that the previous "
                         f"run crashed\n")
            sp.run(f"{scripts_dir}/launch_concat_chr_regions.sh "
                   f"{run_name} {chrom}",
                   shell=True, stdout=sp.PIPE, check=True)
            continue

        if not check_final_output_samples.check_samples(
                run_name, chrom, printing=False):
            status.write(f"The concat job for chrom {chrom} did not produce "
                         f"a vcf with the correct samples in it. "
                         f"Please debug manually\n")
            continue

        # check final number of variants is correct
        job_name = f"check_final_output_num_variants_{chrom}"
        if job_name in job_names:
            status.write(f"Chrom {chrom} is running the final output check "
                         f"num variants job, with name {job_name}\n")
            continue

        error_file = name_from_tag(job_name, error_files)
        if error_file:
            status.write(f"Chrom {chrom} failed the final output num variants"
                         f" check. Errors detailed in file f{error_file}\n")
            continue

        output_file = name_from_tag(f"{job_name}.o", output_files)
        if output_file is None:
            status.write(f"Chrom {chrom} is starting up the final output "
                         f"check num variants job, with name {job_name}\n")
            sp.run(f"{scripts_dir}/launch_check_final_output_num_variants.sh "
                   f"{run_name} {chrom}",
                   shell=True, stdout=sp.PIPE, check=True)
            continue
        else:
            with open(f"{run_dir}/batches/output/{output_file}") as output:
                if "Success" not in output.read():
                    status.write(f"Chrom {chrom} was tested for final output "
                                 f"number of variants but the testing didn't "
                                 f"error or return success. Please "
                                 f"diagnose manually\n")
                    continue

        status.write(f"Chrom {chrom} is finished! Output file "
                     f"{run_dir}/vcfs/chr{chrom}.vcf.gz\n")

def name_from_tag(tag, name_set):
    for name in name_set:
        if tag in name:
            return name
    return None

status_file = (run_dir + "/status_"
               + datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
               + ".txt")
with open(status_file, 'xt', 1) as status:
    chrs_to_check = set(range(1, 23))
    chrs_to_check = check_impute(status,
                                 batch_files,
                                 job_names,
                                 error_files,
                                 output_files,
                                 chrs_to_check)
    chrs_to_check = check_regional_merge(status,
                                         batch_files,
                                         job_names,
                                         error_files,
                                         output_files,
                                         chrs_to_check)
    check_concat(status,
                 final_files,
                 job_names,
                 error_files,
                 output_files,
                 chrs_to_check)

    status.write("""\n\n
How to handle errors:
* To rerun a job which errored, deleted the corresponding file
  (either runs/vcfs/chr<chr>.vcf.gz,
   runs/batches/chr<chr>_pos_<start>_to_<end>.vcf.gz
   or runs/batches/chr<chr>_samples_<start>_to_<end>.vcf.gz)
  and the corresponding .e error file and .o output file
  and rerun this command.
* If a check variants job errored and it seems that the file it
  was checking is corrupted, simply delete that file and the .e
  error file and .o output file on the check variants job
  and rerun this command.
* If a check variants job errored and it seemed like an issue
  with that job and not the underlying file, simply delete the
  .e error file and .o output file for the check variants job
  and rerun this command.

How to handle jobs:
* To get a list of job IDs that are currently running, use qstat -u $(whoami)
* To get a job's full name from it's prefix, run
  qstat -f $(qselect -u $(whoami) | cut -f1 -d. ) | grep <job_prefix>
* To get a job's ID from its name run this command with the job name inserted
  qstat -f $(qselect -u $(whoami) | cut -f1 -d'.') |
                 grep <job_name> -B 1 | head -n 1 | cut -f1 -d.
* To check if a job is hanging, run qpeek <job_id> or qpeek -e <job_id>
  to check what has been written to stdout/err
* To delete a hanging job, run
  qdel <job_id>
""")

print(f"Done. Please see {status_file} for output")
