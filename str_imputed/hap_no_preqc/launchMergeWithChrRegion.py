import subprocess
import sys

if len(sys.argv) != 2:
        print("Expecting exactly 1 argument, the chromosome number")
        exit(-1)
chr = sys.argv[1]

for i in range(1, 250000000, 10000000):
	subprocess.run('qsub -v "INPUT1={},INPUT2={}" mergeWithinChrRegion.pbs'.format(i, chr), shell=True)

