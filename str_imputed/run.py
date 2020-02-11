import subprocess as sp
import argparse
import os, os.path
import glob
import datetime
import sys

sys.path.insert(0, os.environ['UKB'] + "/str_imputed/scripts")

import launch_impute
import check_beagle_output_samples
import check_final_output_samples

parser = argparse.ArgumentParser(description="Continues the imputation process by assessing which jobs need to be run and launching them. Information about the current status of the imputation, issues that need to be resolved, and next steps, to $UKB/str_imputed/runs/<run_name>/status_<date>.txt")
parser.add_argument("run_name", help="The name of this run. Output files will be put in $UKB/str_imputed/runs/<run_name>")
parser.add_argument("--pfile_directory", help="the directory containing the pfiles of the dataset to be imputed (the pfiles must be named chr1 ... chr22). Only necessary the first time running this command")
parser.add_argument("--sample_file", help="the sample file with the list of samples, see $UKB/microarray/*.sample for an example. Only necessary the first time running this command")
parser.add_argument("--readme", help="the description of this run (e.g. what filters were use to create the input .sample file and pfiles, etc.) Only necessary the first time running this command")

args = parser.parse_args()
run_name = args.run_name
pfile_directory = args.pfile_directory
sample_file = args.sample_file
readme = args.readme

def error(msg):
	print(msg, file=sys.stderr)
	exit(-1)

if not TMPDIR in os.environ:
	error("Must set the TMPDIR environment variable")

run_dir = os.environ['UKB'] + "/str_imputed/runs/" + run_name
os.makedirs(run_dir, exist_ok = True)
os.makedirs(f"{run_dir}/batches", exist_ok = True)
os.makedirs(f"{run_dir}/batches/output", exist_ok = True)
os.makedirs(f"{run_dir}/batches/old", exist_ok = True)
os.makedirs(f"{run_dir}/vcfs", exist_ok = True)

info_file = f'{run_dir}/info.txt'
if not os.path.exists(info_file):
	if not pfile_directory or not sample_file:
		error("Options --pfile_directory and --sample_directory are mandatory on first run")


	tmpdir = os.environ['TMPDIR']
	os.path.makedirs(f'{run_dir}', exit_ok=True)
	sp.run(f"plink2 --pfile {pfile_directory}/chr1 --keep {sample_file} --write-samples --out {tmpdir}/{run_name}", stdout = sp.PIPE, shell = True)
	stdout = sp.run(f"wc -l {tmpdir}/{run_name}.id", shell = True, stdout = sp.PIPE, stderr = sp.PIPE).stdout
	num_samples = int(stdout.decode().split()[0]) - 1
	if num_samples < 50000:
		error(f"Error: found too few samples ({num_samples}), something went wrong during pfile loading")	
	print('num_samples', num_samples)
	exit()
	
	with open(info_file, 'w') as info:
		info.write(f"pfile_directory:{pfile_directory}\nsample_file:{sample_file}\nnum_samples:{num_samples}\n")
else:
	with open(info_file) as info:
		fields = {}
		for line in info:
			key, value = line.split(':')
			fields[key.strip()] = value.strip()
		pfile_directory = fields['pfile_directory']
		sample_file = fields['sample_file']
		num_samples = fields['num_samples']

scripts_dir = os.environ['UKB'] + "/str_imputed/scripts"
status_file = run_dir + "/status_" + datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S") + ".txt"

output_files = os.listdir(run_dir + '/batches/output')
errors_output = sp.run(f'{scripts_dir}/check_for_pbs_errors.sh {run_name} | grep "Error in file"', shell=True, universal_newlines=True, stdout=sp.PIPE)
error_files = set()
for line in errors_output.stdout.split('\n'):
	if line.strip() == "":
		continue
	error_files.add(line.split(' ')[3])

job_names = set()
job_ids = sp.run('qselect -u $(whoami)', shell=True, universal_newlines = True, stdout=sp.PIPE).stdout.split('\n')
job_ids = " ".join([job_id.split('.')[0] for job_id in job_ids])
qstat_output = sp.run(f'qstat -f {job_ids} | grep Job_Name ', shell=True, universal_newlines=True, stdout=sp.PIPE).stdout
for line in qstat_output.split("\n"):
	if line.strip() == "":
		continue
	job_names.add(line.split("=")[1].strip())

def check_impute(status, batch_files, job_names, error_files, output_files, chrs_to_check):
	status.write("---Step 1: Imputation---\n")
	keep_checking = set()
	done_with_impute_file = f"{run_dir}/done_with_1_impute.txt"
	with open(done_with_impute_file) as done_with_impute:
		keep_checking.union(int(chrom) for chrom in done_with_impute.read().split('\n'))
	for chrom in chrs_to_check:
		if chrom in keep_checking:
			continue

		#check impute variants
		job_name = f"check_beagle_output_variants_{chrom}"
		if job_name in job_names:
			status.write(f"Chrom {chrom} is running the beagle output check variants job, with name {job_name}.\n")
			check_finished.add(chrom)
			continue

		error_file = name_from_prefix(job_name, error_files)
		if error_file:
			check_finished.add(chrom)
			status.write(f"Chrom {chrom} failed the beagle output variants check. Errors detailed in file f{error_file}.\n")
			continue

		#run launch_impute which handles errors itself
		currently_running_jobs=set()
		for sample_idx in range(1, num_samples + 1, 1000):
			job_name = "impute_STRs_{chrom}_{sample_idx}"
			if job_name in job_names:
				currently_running_jobs.add(job_name)
		ret_code = launch_impute.do_impute(run_name, pfile_directory, sample_file, 
			chrom, readme, command_line=False, ignore_jobs=currently_running_jobs)
		if type(ret_code) == str:
			status.write(f"Launching jobs imputation jobs for chromosome {chrom} resulted in the error {ret_code}\n")
			continue
	
		if ret_code == False:
			status.write(f"Chromosome {chrom} is running imputations")
			continue

		ret_code = check_beagle_output_samples.do_check(run_name, sample_file, chrom, command_line = False)
		if type(ret_code) == str:
			status.write(f"Chromosome {chrom} was imputed with the wrong samples! Please fix this manually. Specific error: {ret_code}\n")
			continue

		output_name = f"check_beagle_output_variants_{chrom}"
		output_file = name_from_prefix(output_name, output_files)
		if output_file is None:
			sp.run(f'{scripts_dir}/check_beagle_output_variants.sh {run_name} {chrom}', shell=True, stdout=sp.PIPE)
		else:
			with open(f"{run_dir}/batches/output/{output_file}") as output:
				if "Success" not in output.read():
					status.write("Chrom {chrom} was tested for imputation output variants but the testing didn't error or return success. Please diagnose manually\n")
					continue

			keep_checking.append(chrom)
		
	with open(done_with_impute_file, 'w'):
		first = True 
		for chrom in keep_checking:
			if not first:
				keep_checking.write("\n")
			first = False
			keep_checking.write(chrom)

	return keep_checking 

'''
def launch_concat(status, chrom):
	status.write(f'Would run {scripts_dir}/launch_concat_chr_regions.sh {run_name} {chrom}\n')
	#sp.run(f'{scripts_dir}/launch_concat_chr_regions.sh {run_name} {chrom}', shell = True, stdout=sp.PIPE)

def launch_check_regional_merge_num_variants(status, chrom, regions = None):
	command = f'Would run {scripts_dir}/launch_concat_chr_regions.sh {run_name} {chrom}'
	if regions is not None:
		command += " --regions "
		for region in regions:
			command += f"{region} "
	status.write(f'Would run command: {command}\n')
	#sp.run(command, shell = True, stdout=sp.PIPE)

def check_regional_merge_num_variants(status, chrom, job_names, error_files, output_files):
	running_notification = False
	error_notification = False
	job_prefix = f"check_merge_output_num_variants_{chrom}"
	regions = []
	for pos in range(1, 5000000, 250000000):
		job_name = f"{job_prefix}_{pos}"
		if job_name in job_names and not running_notification:
			running_notification = True
			status.write(f"The check regional merge num variants jobs chromosome {chrom} are still running. They have names beginning with {job_prefix}\n")

		error_file = name_from_prefix(job_name, error_files)
		if error_file is not None and not error_notification:
			error_notification = True
			status.write(f"Some of the check regional merge num variants jobs for chromosome {chrom} have failed. Please look for error files containing {job_prefix}, such as {error_file}\n")

		output_files = glob.glob(f"{run_dir}/batches/output/{job_name}.o")
		if len(output_files) == 0:
			regions.append(pos)
			continue

		success = False
		for output_file in output_files:
			with open(output_file, 'r') as output:
				if 'Success' in output.read():
					success = True
					break
		if not success:
			regions.append(pos)

	if not running_notification and not error_notification and len(regions) == 0:
		

def check_concat(status, final_files, job_names, error_files, output_files, chrs_to_check):
	status.write("---Step 3: Concat---\n")
	check_finished = set()
	for chrom in chrs_to_check:
		check_regional_merge_num_variants(status, chrom, job_names, error_files, output_files)

		job_name = f"concat_chr_regions_{chrom}"
		if job_name in job_names:
			status.write(f"The concat job for chrom {chrom} is currently running with name {job_name}\n")
			check_regional_merge_num_variants(status, chrom)
			check_finished.add(chrom)
			continue

		error_file = name_from_prefix(job_name, error_files)
		if error_file:
			status.write(f"Concat for chrom {chrom} failed. Errors detailed in file f{error_file}.\n")
			check_regional_merge_num_variants(status, chrom)
			check_finished.add(chrom)
			continue

		tbi_file = f'chr{chrom}.vcf.gz.tbi'
		if not f'chr{chrom}.vcf.gz' in final_files or \
			not tbi_file in final_files:
			continue
			

		if os.stat(run_dir + "/vcfs/" + tbi_file).st_size < 2e4:
			status.write(f"Relaunching the concat job for chrom {chrom} - the tabix index file produced is smaller than anticipated, inticating that the previous run crashed\n")
			launch_concat(chrom)
			check_regional_merge_num_variants(status, chrom)
			check_finished.add(chrom)
			continue

		if not check_final_output_samples.check_samples(run_name, chrom, printing=False):
			status.write(f"Relaunching the concat job for chrom {chrom} - the output vcf file did not have the correct samples, indicating that the previous run failed unexpectedly\n")
			check_regional_merge_num_variants(status, chrom)
			check_finished.add(chrom)
			launch_concat(chrom)
			continue

		#check for length of chrom
		status.write(f"Concat from chrom {chrom} looks good, but need to implement num_variants check\n")
		check_finished.add(chrom)
		
	return check_finished 

def check_regional_merge(status, batch_files, job_names, error_files, output_files, chrs_to_check):
	status.write("---Step 2: Regional Merge---\n")
	check_finished = set()
	for chrom in chrs_to_check:
		ready_for_concat = True
		currently_running = False
		for pos in range(1,250000000,5000000):
			job_name = f"merge_within_region_{chrom}_{pos}"
			end_pos = region+4999999
			if job_name in job_names:
				if not currently_running:
					currently_running = True
					status.write(f"Regional merge jobs for chrom {chrom} are currently running, with names beginning with merge_within_region_{chrom}.\n")
					check_finished.add(chrom)
					ready_for_concat = False
				continue

			error_file = name_from_prefix(job_name, error_files)
			if error_file:
				check_finished.add(chrom)
				ready_for_concat = False
				status.write(f"Regional merge for chrom {chrom} failed for position {pos}-{end_pos}. Errors detailed in file f{error_file}.")

			tbi_file = f'chr{chrom}.vcf.gz.tbi'
			if not f'chr{chrom}.vcf.gz' in final_files or \
				not tbi_file in final_files:
				continue

			if os.stat(run_dir + "/vcfs/" + tbi_file).st_size < 2e5:
				#Simple check that tabix finished
				continue

			if not check_final_output_samples.check_samples(run_name, chrom):
				continue

		#check for length of chrom
		status.write(f"Concat from chrom {chrom} looks good, but need to implement num_variants check\n")
		check_finished.add(chrom)
		
	return check_finished 
'''

def name_from_prefix(prefix, name_set):
	for name in name_set:
		if name.startswith(prefix):
			return name
	return None

with open(status_file, 'xt') as status:
	chrs_to_check = set(range(1,23))
	chrs_to_check = check_impute(status, batch_files, job_names, error_files, output_files, chrs_to_check)

	exit()
	chrs_to_check = list_minus(chrs_to_check, 
		check_regional_merge(status, batch_files, job_names, error_files, output_files, chrs_to_check))
	
	check_impute(status, batch_files, job_names, error_files, output_files, chrs_to_check)

	status.write("""\n\n
How to handle errors:
* To rerun a job which errored, deleted the corresponding file
  (either runs/vcfs/chr<chr>.vcf.gz,
   runs/batches/chr<chr>_pos_<start>_to_<end>.vcf.gz 
   or runs/batches/chr<chr>_samples_<start>_to_<end>.vcf.gz)
  and rerun this command.
* If the error has been resolved and that job ran successfully,
  simply delete the error file and rerun this command.

How to handle jobs:
* To get a list of job IDs that are currently running, use qstat -u $(whoami)
* To get a job's full name from it's prefix, run
  qstat -f $(qselect -u $(whoami) | cut -f1 -d. ) | grep <job_prefix>
* To get a job's ID from its name run this command with the job name inserted
  qstat -f $(qselect -u $(whoami) | cut -f1 -d'.') | grep <job_name> -B 1 | head -n 1 | cut -f1 -d.
* To check if a job is hanging, run qpeek <job_id> or qpeek -e <job_id>
  to check what has been written to stdout/err
* To delete a hanging job, run
  qdel <job_id>
""")

print(f"Please see {status_file} for output")
