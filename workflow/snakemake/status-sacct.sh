#!/usr/bin/env bash

# Check status of Slurm job
# from https://github.com/jdblischak/smk-simple-slurm/blob/main/extras/status-sacct.sh
# necessary for handling silent failures due to timeout/out of mem

jobid="$1"

if [[ "$jobid" == Submitted ]]
then
  echo smk-simple-slurm: Invalid job ID: "$jobid" >&2
  echo smk-simple-slurm: Did you remember to add the flag --parsable to your sbatch call? >&2
  exit 1
fi

## Dont run between 315 and 400 in morning
currenttime=$(date +%H:%M)
if [[ "$currenttime" > "03:15" ]] && [[ "$currenttime" < "04:00" ]]; then
  while [[ $(date +%H:%M) < "04:00" ]]; do sleep 1; done
fi
## End changes

output=`sacct -j "$jobid" --format State --noheader | head -n 1 | awk '{print $1}'`

if [[ $output =~ ^(COMPLETED).* ]]
then
  echo success
elif [[ $output =~ ^(RUNNING|PENDING|COMPLETING|CONFIGURING|SUSPENDED).* ]]
then
  echo running
else
  echo failed
fi
