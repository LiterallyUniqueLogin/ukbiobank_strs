import os, os.path
import sys
import subprocess as sp
import re
from datetime import datetime

if len(sys.argv) != 2:
	print("Expecting exactly 1 argument, the chromosome number")
	exit(-1)
chr = sys.argv[1]

now = datetime.now().strftime("%y_%m_%d_%H_%M_%S")
numSamples = 487409
jobsToRun = [] #list of pairs of (minId, maxId) (inclusive)
jobIds = set()
jobsAlreadyRun = []

def batchName(minId, maxId):
	return "chr{}_samples_{}_to_{}".format(chr, minId, maxId)

def outputDir():
	return "{}/str_imputed/hap_no_preqc/vcf_batches".format(os.environ['UKB'])

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

print("Jobs to run", jobsToRun)
#print("Job ids", jobIds)

existingErrors = sp.run("grep -b10 -v INPUT {}/str_imputed/hap_no_preqc/vcf_batches/output/*.e*".format(os.environ['UKB']),
	shell = True,
	stdout = sp.PIPE)
existingErrors = existingErrors.stdout.decode()

for jobId in [int(match.group(1)) for match in re.finditer("INPUT1 ([0-9]+)", existingErrors)]:
	if jobId not in jobIds:
		print("There's an existing error with job {} but we're not rerunning it. Please solve this problem".format(jobId))
		exit(-1)


existingErrorFiles = sp.run("grep -l -v INPUT {}/str_imputed/hap_no_preqc/vcf_batches/output/*.e*".format(os.environ['UKB']), shell = True, stdout = sp.PIPE)
existingErrorFiles = existingErrorFiles.stdout.decode()
for file in existingErrorFiles.split():
	os.rename(file.replace(".e", ".o"), outputDir() + "/old/" + file.split("/")[-1].replace(".e", ".o"))
	os.rename(file, outputDir() + "/old/" + file.split("/")[-1])

for job in jobsToRun:
	print("Launching job", job)
	sp.run('qsub -v "INPUT1={},INPUT2={},INPUT3={}" impute.pbs'.format(job[0], job[1], chr),
		 shell = True, stdout = sp.PIPE, stderr = sp.PIPE)

