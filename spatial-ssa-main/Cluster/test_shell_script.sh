#!/bin/bash

#SBATCH --job-name="hello_world"
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=1GB

chmod u+x ./test
srun ./test