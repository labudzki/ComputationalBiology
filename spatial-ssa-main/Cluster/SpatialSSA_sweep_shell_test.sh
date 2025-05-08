#!/bin/bash

#SBATCH --job-name="Yeast_polarization_sweep_phase_seperation"
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=1GB
#SBATCH --output=./Data/Output/O_%x_%A_%a.out
#SBATCH --array=1-125

# Read the ith line of the sweep file and put into sweep setup
# https://stackoverflow.com/questions/6022384/bash-tool-to-get-nth-line-from-a-file
# https://stackoverflow.com/questions/6749128/store-output-of-sed-into-a-variable
sweep_file="sweep_setups_ps.txt"
sweep_setup=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sweep_file})

echo "Task ID"
echo "${SLURM_ARRAY_TASK_ID}"

echo "Array job ID"
echo "${SLURM_ARRAY_JOB_ID}"

echo "Array task ID"
echo "${SLURM_ARRAY_TASK_ID}"


echo "Task setup"
echo "${sweep_setup}"

timestamp=$(date +%Y-%m-%d_%H-%M-%S)
dirname="Data/${timestamp}_job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Making directory ${dirname}"
mkdir -p "${dirname}"

# chmod u+x ./SpatialSSA_setup_ps
# srun ./SpatialSSA_setup_ps ${dirname} ${sweep_setup}