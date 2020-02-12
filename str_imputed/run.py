import subprocess as sp
import argparse
import os, os.path
import glob
import datetime
import sys

sys.path.insert(0, os.environ['UKB'] + "/str_imputed/scripts")

import launch_impute
import check_beagle_output_samples
import launch_merge_within_region
import check_merge_output_samples
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

if 'TMPDIR' not in os.environ:
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
	os.makedirs(f'{run_dir}', exist_ok=True)
	sp.run(f"plink2 --pfile {pfile_directory}/chr1 --keep {sample_file} --write-samples --out {tmpdir}/{run_name}", stdout = sp.PIPE, shell = True)
	stdout = sp.run(f"wc -l {tmpdir}/{run_name}.id", shell = True, stdout = sp.PIPE, stderr = sp.PIPE).stdout
	num_samples = int(stdout.decode().split()[0]) - 1
	if num_samples < 50000:
		error(f"Error: found too few samples ({num_samples}), something went wrong during pfile loading")	
	
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
		num_samples = int(fields['num_samples'])

scripts_dir = os.environ['UKB'] + "/str_imputed/scripts"
final_files =  os.listdir(run_dir + '/vcfs')
batch_files = os.listdir(run_dir + '/batches')
output_files = os.listdir(run_dir + '/batches/output')

print("Collecting errors")
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
	if os.path.exists(done_with_impute_file):
		with open(done_with_impute_file) as done_with_impute:
			#TODO fix this
			keep_checking.union(int(chrom) for chrom in done_with_impute.read().split('\n') if chrom != "")

	for chrom in chrs_to_check:
		if chrom in keep_checking:
			continue

		print(f"Checking chrom {chrom} impute progress")	

		#run launch_impute which handles errors itself
		currently_running_jobs=set()
		for sample_idx in range(1, num_samples + 1, 1000):
			job_name = "impute_STRs_{chrom}_{sample_idx}"
			if job_name in job_names:
				currently_running_jobs.add(job_name)
		ret_code = launch_impute.do_impute(run_name, pfile_directory, sample_file, 
			chrom, readme, command_line=False, ignore_jobs=currently_running_jobs)
		if type(ret_code) == str:
			status.write(f"Launching imputation jobs for chromosome {chrom} resulted in the error {ret_code}\n")
			continue
	
		if ret_code == False:
			status.write(f"Chromosome {chrom} is running imputations with names beginning with impute_STRs_{chrom}\n")
			continue

		ret_code = check_beagle_output_samples.do_check(run_name, sample_file, chrom, command_line = False)
		if type(ret_code) == str:
			status.write(f"Chromosome {chrom} was imputed with the wrong samples! Please fix this manually. Specific error: {ret_code}\n")
			continue
		
		#check imputation variants are correct
		job_name = f"check_beagle_output_variants_{chrom}"
		if job_name in job_names:
			status.write(f"Chrom {chrom} is running the beagle output check variants job, with name {job_name}\n")
			continue

		error_file = name_from_prefix(job_name, error_files)
		if error_file:
			status.write(f"Chrom {chrom} failed the beagle output variants check. Errors detailed in file f{error_file}\n")
			continue

		output_file = name_from_prefix(f"{job_name}.o", output_files)
		if output_file is None:
			status.write(f"Chrom {chrom} is starting up the beagle output check variants job, with name {job_name}\n")
			sp.run(f'{scripts_dir}/launch_check_beagle_output_variants.sh {run_name} {chrom}', shell=True, stdout=sp.PIPE)
		else:
			output_file = f"{run_dir}/batches/output/{output_file}"
			with open(output_file) as output:
				if "Success" not in output.read():
					status.write(f"Chrom {chrom} was tested for imputation output variants but the testing didn't error or return success. Please diagnose manually. See file {output_file}\n")
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

def check_regional_merge(status, batch_files, job_names, error_files, output_files, chrs_to_check):
	status.write("---Step 2: Regional Merge---\n")
	keep_checking = set()
	done_with_regional_merge_file = f"{run_dir}/done_with_2_regional_merge.txt"
	if os.path.exists(done_with_regional_merge_file):
		with open(done_with_regional_merge_file) as done_with_regional_merge:
			keep_checking.union(int(chrom) for chrom in done_with_regional_merge.read().split('\n') if chrom != "")

	for chrom in chrs_to_check:
		if chrom in keep_checking:
			continue
		
		print(f"Checking chrom {chrom} regional merge progress")	

		keep_checking_this_chrom = True

		#run launch_impute which handles errors itself
		sample_failures = check_merge_output_samples.do_check(run_name, chrom, command_line = False)
		currently_running_jobs=False
		errored_jobs=False
		final_sample_failures = set()
		new_jobs=set()
		check_num_variants=False
		check_num_variants_errors=False
		new_checks=set()
		for pos in range(1, 250000000, 5000000):
			job_name = "merge_within_region_{chrom}_{pos}"
			if job_name in job_names:
				currently_running_jobs = True
				keep_checking_this_chrom = False
				continue
			temp_error_file = name_from_prefix(job_name, error_files)
			if temp_error_file is not None:
				error_file = temp_error_file
				errored_jobs = True
				keep_checking_this_chrom = False
				continue

			batch_file = name_from_prefix("chr{chrom}_pos_{pos}", batch_files)
			if batch_file is None:
				new_jobs.add(pos)
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
			
			temp_check_error_file = name_from_prefix(job_name, error_files)
			if temp_check_error_file is not None:
				check_error_file = temp_check_error_file
				keep_checking_this_chrom = False
				check_num_variants_errors= True
				continue

			output_file = name_from_prefix(f"{job_name}.o", output_files)
			if output_file is None:
				new_checks.add(pos)
				keep_checking_this_chrom = False
			else:
				with open(f"{run_dir}/batches/output/{output_file}") as output:
					if 'Success' not in output.read():
						status.write(f"Chromosome {chrom} check regional merge num varaints jobs finished without either error or success. Please debug manually. See file {output_file}\n")
						keep_check_this_chrom = False
				#this position is fully successful, don't need to do anything
						

		if currently_running_jobs:
			status.write(f"Chromosome {chrom} is running regional merges\n")
		if errored_jobs:
			status.write(f"Chromosome {chrom} failed the regional merge in at least one region. Errors detailed in files named similarly to f{error_file}\n")

		if len(new_jobs) > 0:
			regions = " ".join(list(new_jobs))
			msg = launch_merge_within_region.do_launch(run_name, chrom, regions, command_line=False)
			if msg is not None:
				status.write(f"Chromosome {chrom} failed to launch regional merge jobs. Please debug manually. See error {msg}\n")
			else:
				status.write(f"Chromsome {chrom} is launching new regional merge jobs\n")
		
		if len(final_sample_failures) != 0:
			status.write(f"Chromosome {chrom} was merged with the wrong samples in regions {final_sample_failures}! Please fix this manually. Specific error: {ret_code}\n")
	
		if check_num_variants:
			status.write(f"Chrom {chrom} is running merge output check variants job(s), with names beginning with check_merge_output_num_varaints_{chrom}\n")
		if check_num_variants_errors:
			status.write(f"Chromosome {chrom} failed the regional merge num variants check in at least one region. Errors detailed in files named similarly to f{check_error_file}\n")

		if len(new_checks) > 0:
			regions = " ".join(list(new_checks))
			sp.run(f'{scripts_dir}/launch_check_merge_output_num_variants.sh {run_name} {chrom} {regions}', shell=True, stdout=sp.PIPE)
			status.write(f"Chromosome {chrom} is running new regional merge num variants checks. Job names similar to check_merge_output_num_variants_{chrom}\n")

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

def check_concat(status, final_files, job_names, error_files, output_files, chrs_to_check):
	status.write("---Step 3: Concat---\n")
	for chrom in chrs_to_check:

		print(f"Checking chrom {chrom} concat progress")	

		job_name = f"concat_chr_regions_{chrom}"
		if job_name in job_names:
			status.write(f"The concat job for chrom {chrom} is currently running with name {job_name}\n")
			continue

		error_file = name_from_prefix(job_name, error_files)
		if error_file:
			status.write(f"Concat for chrom {chrom} failed. Errors detailed in file f{error_file}\n")
			continue

		tbi_file = f'chr{chrom}.vcf.gz.tbi'
		if not f'chr{chrom}.vcf.gz' in final_files or \
			not tbi_file in final_files:
			status.write(f"Launching the concat job for chrom {chrom}, job name concat_chr_regions_{chrom}\n")
			sp.run(f'{scripts_dir}/launch_concat_chr_regions.sh {run_name} {chrom}', shell = True, stdout=sp.PIPE)
			continue

		if os.stat(run_dir + "/vcfs/" + tbi_file).st_size < 2e4:
			status.write(f"Relaunching the concat job for chrom {chrom} - the tabix index file produced is smaller than anticipated, inticating that the previous run crashed\n")
			sp.run(f'{scripts_dir}/launch_concat_chr_regions.sh {run_name} {chrom}', shell = True, stdout=sp.PIPE)
			continue

		if not check_final_output_samples.check_samples(run_name, chrom, printing=False):
			status.write(f"The concat job for chrom {chrom} did not produce a vcf with the correct samples in it. Please debug manually\n")
			continue

		#check final number of variants is correct
		job_name = f"check_final_output_num_variants_{chrom}"
		if job_name in job_names:
			status.write(f"Chrom {chrom} is running the final output check num variants job, with name {job_name}\n")
			continue

		error_file = name_from_prefix(job_name, error_files)
		if error_file:
			status.write(f"Chrom {chrom} failed the final output num variants check. Errors detailed in file f{error_file}\n")
			continue

		output_file = name_from_prefix(f"{job_name}.o", output_files)
		if output_file is None:
			status.write(f"Chrom {chrom} is starting up the final output check num variants job, with name {job_name}\n")
			sp.run(f'{scripts_dir}/launch_check_final_output_num_variants.sh {run_name} {chrom}', shell=True, stdout=sp.PIPE)
			continue
		else:
			with open(f"{run_dir}/batches/output/{output_file}") as output:
				if "Success" not in output.read():
					status.write(f"Chrom {chrom} was tested for final output number of variants but the testing didn't error or return success. Please diagnose manually\n")
					continue
		
		status.write(f"Chrom {chrom} is finished! Output file {run_dir}/vcfs/chr{chrom}.vcf.gz\n")

def name_from_prefix(prefix, name_set):
	for name in name_set:
		if name.startswith(prefix):
			return name
	return None

status_file = run_dir + "/status_" + datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S") + ".txt"
with open(status_file, 'xt', 1) as status:
	chrs_to_check = set(range(1,23))
	chrs_to_check = check_impute(status, batch_files, job_names, error_files, output_files, chrs_to_check)
	chrs_to_check = check_regional_merge(status, batch_files, job_names, error_files, output_files, chrs_to_check)
	check_concat(status, final_files, job_names, error_files, output_files, chrs_to_check)

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

print(f"Done. Please see {status_file} for output")
