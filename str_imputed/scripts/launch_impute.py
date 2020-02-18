import argparse
import glob
import os
import os.path
import re
import subprocess as sp
import sys
from datetime import datetime

UKB = os.environ['UKB']
COMMAND_LINE = False

def error(msg):
    if COMMAND_LINE:
        print(msg, file=sys.stderr)
        exit(-1)
    else:
        return msg

def set_up_env(run_name, chrom, readme):
    os.makedirs(f"{UKB}/str_imputed/runs/{run_name}/", exist_ok=True)
    os.makedirs(f"{UKB}/str_imputed/runs/{run_name}/batches", exist_ok=True)
    os.makedirs(f"{UKB}/str_imputed/runs/{run_name}/batches/output",
                exist_ok=True)
    os.makedirs(f"{UKB}/str_imputed/runs/{run_name}/batches/old",
                exist_ok=True)

    # Write the README for this run, or check that it exists as specified
    if not os.path.exists(f"{UKB}/str_imputed/runs/{run_name}/README"):
        if readme is None or readme == "":
            error("Error: No README for this run already exists, but didn't "
                  "specify one with the --readme argument.")
        else:
            with open(f"{UKB}/str_imputed/runs/{run_name}/README",
                      'w') as readme_file:
                readme_file.write(readme + "\n")
            with open(f"{UKB}/str_imputed/run_readmes/{run_name}_README",
                      'w') as readme_file:
                readme_file.write(readme + "\n")
    else:
        if readme is not None and readme != "":
            with open(f"{UKB}/str_imputed/runs/{run_name}/README",
                      'r') as readme_file:
                if (readme + "\n") != readme_file.read():
                    error("Error: Found a different description in the README"
                          " file than the currently intended one. Either "
                          "delete the README file so a new one can be "
                          "written, or remove the --readme flag")


# if COMMAND_LINE==False
# this method returns True if the chromosome is done with imputation
# (no jobs were launched) and False if jobs were successfully launched
# and returns a string error message if there was a problem
def do_impute(run_name,
              pfile_dir,
              sample_file,
              chrom,
              readme,
              ignore_jobs=set()):
    pfile_no_ext = f"{pfile_dir}/chr{chrom}"
    if not os.path.exists(f"{pfile_no_ext}.pgen"):
        error(f"Error: expected file {pfile_no_ext}.pgen to exist")

    if "TMPDIR" not in os.environ:
        error("Error, expected the TMPDIR environment variable to be set.")

    set_up_env(run_name, chrom, readme)

    tmpdir = os.environ['TMPDIR']

    now = datetime.now().strftime("%y_%m_%d_%H_%M_%S")

    # get the total number of samples
    sp.run(f"plink2 --pfile {pfile_no_ext} --keep {sample_file} "
           f"--write-samples --out {tmpdir}/{run_name}",
           stdout=sp.PIPE, shell=True, check=True,)
    stdout = sp.run(f"wc -l {tmpdir}/{run_name}.id",
        shell=True, stdout=sp.PIPE, stderr=sp.PIPE, check=True).stdout
    numSamples = int(stdout.decode().split()[0]) - 1
    if numSamples < 50000:
        error(f"Error: found too few samples ({numSamples}), something went "
              "wrong during pfile loading")

    # Figure out which jobs have already been run
    jobsToRun = []  # list of pairs of (minId, maxId) (inclusive)
    job_ids = set()

    def batchName(minId, maxId):
        return "chr{}_samples_{}_to_{}".format(chrom, minId, maxId)

    def outputDir():
        return f"{UKB}/str_imputed/runs/{run_name}/batches"

    def outputLocNoExt(minId, maxId):
        return "{}/{}".format(outputDir(), batchName(minId, maxId))

    def outputLoc(minId, maxId):
        return outputLocNoExt(minId, maxId) + ".vcf.gz.tbi"

    def outputLogLoc(minId, maxId):
        return outputLocNoExt(minId, maxId) + ".log"

    for minId in range(1, numSamples, 1000):
        if minId in ignore_jobs:
            continue
        maxId = min(minId + 999, numSamples)

        if not os.path.exists(outputLoc(minId, maxId)):
            jobsToRun.append((minId, maxId))
            job_ids.add(minId)
            continue

        if os.path.getsize(outputLoc(minId, maxId)) < 10000:
            jobsToRun.append((minId, maxId))
            job_ids.add(minId)
            continue

        if os.path.exists(outputLogLoc(minId, maxId)):
            os.rename(outputLogLoc(minId, maxId), outputDir() + "/output/"
                      + batchName(minId, maxId) + ".log")

    for job in jobsToRun:
        minId = job[0]
        maxId = job[1]
        if os.path.exists(outputLocNoExt(minId, maxId) + ".vcf.gz"):
            os.rename(outputLocNoExt(minId, maxId) + ".vcf.gz", outputDir()
                      + "/old/" + now + "_"
                      + batchName(minId, maxId) + ".vcf.gz")
        if os.path.exists(outputLoc(minId, maxId)):
            os.rename(outputLoc(minId, maxId), outputDir()
                      + "/old/" + now + "_"
                      + batchName(minId, maxId) + ".vcf.gz.tbi")
        if os.path.exists(outputLogLoc(minId, maxId)):
            os.rename(outputLogLoc(minId, maxId), outputDir()
                      + "/old/" + now + "_"
                      + batchName(minId, maxId) + ".log")
    if COMMAND_LINE:
        print(f"chr{chrom}: Jobs to run", jobsToRun)

    if len(jobsToRun) == 0:
        if COMMAND_LINE:
            exit()
        else:
            return len(ignore_jobs) == 0

    existingErrorIds = set()
    existingErrorFiles = set()
    for file_name in glob.glob(
            f"{UKB}/str_imputed/runs/{run_name}/batches/output/*.e*"):
        with open(file_name) as file:
            contents = file.readlines()
            error = False
            current_chrom = False
            for line in contents:
                if f"INPUT3 {chrom}" in line:
                    current_chrom = True
                if "INPUT1" not in line:
                    error = True
            if error and current_chrom:
                existingErrorFiles.add(file_name)
                matches = re.findall("INPUT1 ([0-9]+)", " ".join(contents))
                if len(matches) != 1:
                    error("Found more than 1 INPUT1 in an error file, "
                          "confused")
                existingErrorIds.add(int(matches[0]))

    fix_jobs = set()
    for job_id in existingErrorIds:
        if job_id not in job_ids:
            fix_jobs.add(job_id)

    if len(fix_jobs) > 0:
        error(f"There are existing errors with jobs {fix_jobs} "
              f"but we're not rerunning them. Please solve this "
              f"problem by either removing the associated .vcf.gz.tbi "
              f"files from {UKB}/str_imputed/runs/{run_name}/batches, "
              f"which will cause the jobs to be rerun, or remove the "
              f"associated .e error files from "
              f"{UKB}/str_imputed/runs/{run_name}/batches/output, "
              f"which will cause the job not to be rerun "
              f"(the error file can be found by grepping the *.e* "
              f"files in that directory for 'INPUT1 <job_id> '")

    # We're going to rerun the files which have existing errors,
    # so move all the errors to the old directory
    for file in existingErrorFiles:
        os.rename(file.replace(".e", ".o"),
                  outputDir() + "/old/"
                  + file.split("/")[-1].replace(".e", ".o"))
        os.rename(file, outputDir() + "/old/" + file.split("/")[-1])

    for job in jobsToRun:
        # create munged impute.pbs file
        sp.run(f"sed -e 's/%RUN_NAME%/{run_name}/g' "
               "-e 's/%CHROM%/{chrom}/g' "
               "-e 's/%SAMPLE%/{job[0]}/g' "
               "{UKB}/str_imputed/scripts/impute.pbs > {tmpdir}/impute_{run_name}.pbs",
                check=True, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        if COMMAND_LINE:
            print(f"chr{chrom}: Launching job {job}")
        sp.run(f"qsub -v "
               f"INPUT1={job[0]},INPUT2={job[1]},INPUT3={chrom},"
               f"INPUT4={pfile_dir},INPUT5={sample_file}"
               f"{tmpdir}/impute_{run_name}.pbs",
               check=True, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_name",
                        help="Output files will be put in "
                             "$UKB/str_imputed/runs/run_name/batches")
    parser.add_argument("pfile_directory",
                        help="the directory containing the pfiles of the "
                             "dataset to be imputed (the pfiles must be named "
                             "chr1 ... chr22)")
    parser.add_argument("sample_file",
                        help="the .sample file with the list of samples, see "
                             "$UKB/microarray/*.sample for an example")
    parser.add_argument("chromosome_number",
                        help="the number of the chromosome to impute")
    parser.add_argument("--readme",
                        default="",
                        help="the description of this run (e.g. what "
                             "filters were use to create the input ."
                             "sample file and pfiles, etc.) Required if "
                             "a README does not already exist for this run.")

    args = parser.parse_args()

    run_name = args.run_name
    pfile_dir = args.pfile_directory
    sample_file = args.sample_file
    chrom = args.chromosome_number
    readme = args.readme
    do_impute(run_name, pfile_dir, sample_file, chrom, readme)

if __name__ == "__main__":
    COMMAND_LINE = True
    main()
