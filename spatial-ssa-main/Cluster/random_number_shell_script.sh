#!/bin/bash

#SBATCH --job-name="hello_world"
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=1GB

timestamp=$(date +%Y-%m-%d_%H-%M-%S)
dirname="Data/${timestamp}_job_$SLURM_JOB_ID"
echo "Making directory ${dirname}"
mkdir -p "${dirname}"

chmod u+x ./test
srun ./test "${dirname}"