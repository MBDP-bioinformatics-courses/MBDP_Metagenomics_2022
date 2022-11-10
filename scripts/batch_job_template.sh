#!/bin/bash -l
#SBATCH --job-name=XX
#SBATCH --account=project_XXX
#SBATCH --error=00_LOGS/XX-%j.err
#SBATCH --output=00_LOGS/XX-%j.out
#SBATCH --time=XX:XX:XX
#SBATCH --mem=XXG
#SBATCH --partition=small
#SBATCH --cpus-per-task=XX
#SBATCH --nodes=1
#SBATCH --ntasks=1

# load needed module
module load XXX

# use this to specify the number of threads 
$SLURM_CPUS_PER_TASK
