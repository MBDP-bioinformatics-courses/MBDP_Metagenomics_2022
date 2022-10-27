#!/bin/bash -l
#SBATCH --job-name=XX
#SBATCH --account=project_XXX
#SBATCH --error=slurm_logs/XX-%j.err
#SBATCH --output=slurm_logs/XX-%j.out
#SBATCH --time=XX:XX:XX
#SBATCH --mem=XXG
#SBATCH --partition=small
#SBATCH --cpus-per-task=XX
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load XXX

# use this to specify the number of threads 
$SLURM_CPUS_PER_TASK
