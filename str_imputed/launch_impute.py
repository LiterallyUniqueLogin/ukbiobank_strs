import os, os.path
import sys
import subprocess as sp
import re
from datetime import datetime
import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument("run_name", help="Output files will be put in $UKB/str_imputed/run_name")
parser.add_argument("pfile_directory", help="the directory containing the pfiles of the dataset to be imputed (the pfiles must be named chr1 ... chr22)")
parser.add_argument("sample_directory", help="the directory containing the .sample file with the list of samples, see $UKB/microarray/*.sample for an example")
parser.add_argument("chromosome_number", help="the number of the chromosome to impute")
parser.add_argument("--readme", help="the description of this run (e.g. what filters were use to create the input .sample file and pfiles, etc.) Required if a README does not already exist for this run.", default="")

args= parser.parse_args()

run_name = args.run_name
pfile_dir = args.pfile_directory
sample_dir = args.sample_directory
chr = args.chromosome_number
readme = args.readme

pfile_no_ext = f"{pfile_dir}/chr{chr}"
if not os.path.exists(f"{pfile_no_ext}.pgen"):
	print(f"Error: expected file {pfile_no_ext}.pgen to exist", file = sys.stderr)
	exit(-1)
if len(glob.glob(f"{sample_dir}/*.sample")) != 1:
	print(f"Error: expected exactly one file matching {sample_dir}/*.sample", file = sys.stderr)
	exit(-1)

ukb = os.environ['UKB']
if not "TMPDIR" in os.environ:
	print("Error, expected the TMPDIR environment variable to be set.", file = sys.stderr)
	exit(-1)
tmpdir = os.environ['TMPDIR']

os.makedirs(f"{ukb}/str_imputed/{run_name}/", exist_ok = True)
os.makedirs(f"{ukb}/str_imputed/{run_name}/batches", exist_ok = True)
os.makedirs(f"{ukb}/str_imputed/{run_name}/batches/output", exist_ok = True)
os.makedirs(f"{ukb}/str_imputed/{run_name}/batches/old", exist_ok = True)

#Write the README for this run, or check that it exists as specified
if not os.path.exists(f"{ukb}/str_imputed/{run_name}/README"):
	if readme == "":
		print("Error: No README for this run already exists, but didn't specify one with the --readme argument.", file = sys.stderr)
		exit(-1)
	else:
		with open(f"{ukb}/str_imputed/{run_name}/README", 'w') as README_file:
			README_file.write(readme + "\n")
else:
	if readme != "":
		with open(f"{ukb}/str_imputed/{run_name}/README") as README_file:
			if (readme + "\n") != README_file.read():
				print("Error: Found a different description in the README file than the currently intended one. Either delete the README file so a new one can be written, or remove the --readme flag", file = sys.stderr)
				exit(-1)

now = datetime.now().strftime("%y_%m_%d_%H_%M_%S")

#get the total number of samples
sp.run(f"plink2 --pfile {pfile_no_ext} --write-samples --out {tmpdir}/{run_name}",
	stdout = sp.PIPE,
	shell = True)
stdout = sp.run(f"wc -l {tmpdir}/{run_name}.id", 
	shell = True, stdout = sp.PIPE, stderr = sp.PIPE).stdout
numSamples = int(stdout.decode().split()[0]) - 1
if numSamples < 50000:
	print(f"Error: found too few samples ({numSamples}), something went wrong during pfile loading", file = sys.stderr)
	exit(-1)

#Figure out which jobs have already been run
jobsToRun = [] #list of pairs of (minId, maxId) (inclusive)
jobIds = set()
jobsAlreadyRun = []

def batchName(minId, maxId):
	return "chr{}_samples_{}_to_{}".format(chr, minId, maxId)

def outputDir():
	return f"{ukb}/str_imputed/{run_name}/batches"

def outputLocNoExt(minId, maxId):
	return "{}/{}".format(outputDir(), batchName(minId, maxId))

def outputLoc(minId, maxId):
	return outputLocNoExt(minId, maxId) + ".vcf.gz.tbi"

def outputLogLoc(minId, maxId):
	return outputLocNoExt(minId, maxId) + ".log"

for minId in range(1, numSamples, 1000):
	maxId = min(minId+999, numSamples)
	
	if not os.path.exists(outputLoc(minId, maxId)):
		jobsToRun.append((minId, maxId))
		jobIds.add(minId)
		continue

	if os.path.getsize(outputLoc(minId, maxId)) < 10000:
		jobsToRun.append((minId, maxId))
		jobIds.add(minId)
		continue

	jobsAlreadyRun.append((minId, maxId))
	if os.path.exists(outputLogLoc(minId, maxId)):
		os.rename(outputLogLoc(minId, maxId), outputDir() + "/output/" + batchName(minId, maxId) + ".log")

for job in jobsToRun:
	minId = job[0]
	maxId = job[1]
	if os.path.exists(outputLogLoc(minId, maxId)):
		os.rename(outputLogLoc(minId, maxId), outputDir() + "/old/" + now + "_" + batchName(minId, maxId) + ".log")
	if os.path.exists(outputLocNoExt(minId, maxId) + ".vcf.gz"):
		os.rename(outputLocNoExt(minId, maxId) + ".vcf.gz", outputDir() + "/old/" + now + "_"  + batchName(minId, maxId) + ".vcf.gz")
	if os.path.exists(outputLoc(minId, maxId)):
		os.rename(outputLoc(minId, maxId), outputDir() + "/old/" + now + "_"  + batchName(minId, maxId) + ".vcf.gz.tbi")

print(f"chr{chr}: Jobs to run", jobsToRun)

if len(jobsToRun) == 0:
	exit()

existingErrorIds = set()
existingErrorFiles = set()
for file_name in glob.glob(f"{ukb}/str_imputed/{run_name}/batches/output/*.e*"):
	with open(file_name) as file:
		contents = file.readlines()
		error = False
		current_chrom = False
		for line in contents:
			if f"INPUT3 {chr}" in line:
				current_chrom = True
			if not "INPUT1" in line:
				error = True
		if error and current_chrom:
			existingErrorFiles.add(file_name)
			matches = re.findall("INPUT1 ([0-9]+)", " ".join(contents))
			if len(matches) != 1:
				print("Found more than 1 INPUT1 in an error file, confused", file = sys.stderr)
				exit(-1)
			existingErrorIds.add(int(matches[0]))

fix_jobs = set()
for jobId in existingErrorIds:
	if jobId not in jobIds:
		fix_jobs.add(jobId)

if len(fix_jobs) > 0:
	print(f"There are existing errors with jobs {fix_jobs} but we're not rerunning them. Please solve this problem by either removing the associated .vcf.gz.tbi files from {ukb}/str_imputed/{run_name}/batches, which will cause the jobs to be rerun, or remove the associated .e error files from {ukb}/str_imputed/{run_name}/batches/output, which will cause the job not to be rerun (the error file can be found by grepping the *.e* files in that directory for 'INPUT1 <job_id> '", file = sys.stderr)
	exit(-1)

#create munged impute.pbs file
sp.run(f"sed -e 's/%RUN_NAME%/{run_name}/g' {ukb}/str_imputed/impute.pbs > {tmpdir}/impute_{run_name}.pbs",
	shell = True, stdout = sp.PIPE, stderr = sp.PIPE)

#We're going to rerun the files which have existing errors,
#so move all the errors to the old directory
for file in existingErrorFiles:
	os.rename(file.replace(".e", ".o"), outputDir() + "/old/" + file.split("/")[-1].replace(".e", ".o"))
	os.rename(file, outputDir() + "/old/" + file.split("/")[-1])

for job in jobsToRun:
	print(f"chr{chr}: Launching job {job}")
	sp.run(f'qsub -v "INPUT1={job[0]},INPUT2={job[1]},INPUT3={chr},INPUT4={pfile_dir},INPUT5={sample_dir}" {tmpdir}/impute_{run_name}.pbs',
		 shell = True, stdout = sp.PIPE, stderr = sp.PIPE)

