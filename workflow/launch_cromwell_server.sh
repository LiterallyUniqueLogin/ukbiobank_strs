#!/usr/bin/env sh

#SBATCH --job-name=cromwel_server
#SBATCH --account=ddp268
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-oneweek
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=05:23:30
#SBATCH --no-requeue
#SBATCH --export=ALL
#SBATCH --output=cromwell-workflow-logs/server.o%j.%N

declare -xr GALYLEO_LAUNCH_DIR="${SLURM_SUBMIT_DIR}"
declare -xr LOCAL_SCRATCH_DIR="/scratch/${USER}/job_${SLURM_JOB_ID}"

declare -xi PORT=-1
declare -xir LOWEST_EPHEMERAL_PORT=49152
declare -i random_ephemeral_port=-1

module purge
module load cpu

while (( "${PORT}" < 0 )); do
  while (( "${random_ephemeral_port}" < "${LOWEST_EPHEMERAL_PORT}" )); do
    random_ephemeral_port="$(od -An -N 2 -t u2 -v < /dev/urandom)"
  done
  ss -nutlp | cut -d : -f2 | grep "^${random_ephemeral_port})" > /dev/null
  if [[ "${?}" -ne 0 ]]; then
    PORT="${random_ephemeral_port}"
  fi
done

curl -s "https://manage.${GALYLEO_REVERSE_PROXY_FQDN}/linktoken.cgi?token=${REVERSE_PROXY_TOKEN}\&jobid=${SLURM_JOB_ID}"


java -Dconfig.file=workflow/cromwell.conf -Dwebservice.port="${PORT}" -jar utilities/cromwell-86-90af36d-SNAP-600ReadWaitTimeout.jar server
if [[ "${?}" -ne 0 ]]; then
  echo 'ERROR: Failed to launch Cromwell.'
  exit 1
fi

curl "https://manage.${GALYLEO_REVERSE_PROXY_FQDN}/redeemtoken.cgi?token=${REVERSE_PROXY_TOKEN}&port=${PORT}"

wait

curl "https://manage.${GALYLEO_REVERSE_PROXY_FQDN}/destroytoken.cgi?token=${REVERSE_PROXY_TOKEN}"
