import os, os.path
import shutil
import sys
import subprocess as sp

if len(sys.argv) != 2:
	print("Expecting exactly 1 argument, the chromosome number")
	exit(-1)
chr = sys.argv[1]

numSamples = 487409
jobsToRun = set()
jobsAlreadyRun = set()

def batchName(minId, maxId):
	return "chr{}_samples_{}_to_{}".format(chr, minId, maxId)

def outputDir():
	return "{}/str_imputed/hap_no_preqc/vcf_batches".format(os.environ['UKB'])

def outputLocNoExt(minId, maxId):
	return "{}/{}".format(outputDir(), batchName(minId, maxId))

def outputVCF(minId, maxId):
	return outputLocNoExt(minId, maxId) + ".vcf.gz"

def outputTBI(minId, maxId):
	return outputLocNoExt(minId, maxId) + ".vcf.gz.tbi"

#Figure out which jobs to run
#Run jobs in one of the following five scenarios:
#The header hasn't yet been added to the vcf file
#the tbi file doesn't exist for the new vcf file
#the tbi file exists, but is older than the new vcf file
#the tbi file exists but is too small
#(either of the last two imply the program crashed in the middle of creating the new vcf file 
# or new tbi file and never)
for minId in range(1, numSamples, 1000):
	maxId = min(minId+999, numSamples)

	if not os.path.exists(outputVCF(minId, maxId)):
		shutil.copyfile(outputDir() + "/cache/" + batchName(minId, maxId) + ".vcf.gz",
				outputVCF(minId, maxId))

	output = sp.run("zcat {} | head -n 30 | grep AP1".format(outputVCF(minId, maxId)),
		 shell = True, stdout = sp.PIPE, stderr = sp.PIPE)

	if not output:
		jobsToRun.add(minId)
		continue

	if not os.path.exists(outputTBI(minId, maxId)):
		jobsToRun.add(minId)
		continue

	if os.path.getmtime(outputTBI(minId, maxId)) < os.path.getmtime(outputVCF(minId, maxId)):
		jobsToRun.add(minId)
		continue

	if os.path.getsize(outputTBI(minId, maxId)) < 10000:
		jobsToRun.add(minId)
		continue

print("Jobs to run", jobsToRun)

for job in jobsToRun:
	print("Launching job", job)
	sp.run('qsub -v "INPUT1={},INPUT2={}" add_ap_header.pbs'.format(job, chr),
		 shell = True, stdout = sp.PIPE, stderr = sp.PIPE)

