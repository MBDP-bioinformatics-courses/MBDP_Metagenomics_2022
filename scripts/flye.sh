#!/bin/bash -l
#SBATCH --job-name assembly
#SBATCH --output 00_LOGS/assembly-%j.out
#SBATCH --error 00_LOGS/assembly-%j.err
#SBATCH --time 24:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 40
#SBATCH --mem 300G
#SBATCH --account project_2001499
#SBATCH --gres=nvme:200

/scratch/project_2001499/envs/Flye/bin/flye  \
        --nano-raw 02_TRIMMED_DATA/SRR11673980_chop.fastq.gz \
        --threads $SLURM_CPUS_PER_TASK \
        --meta \
        --out-dir 03_ASSEMBLY