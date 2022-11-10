#!/bin/bash -l
#SBATCH --job-name=polish
#SBATCH --account=project_2001499
#SBATCH --error=00_LOGS/polish-%j.err
#SBATCH --output=00_LOGS/polish-%j.out
#SBATCH --time=2-00:00:00
#SBATCH --mem=150GG
#SBATCH --partition=small
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load bowtie2/2.4.4 
bowtie2-build 03_ASSEMBLY/assembly.fasta \
                04_POLISH/assembly

bowtie2 -1 02_TRIIMMED_DATA/SRR11674041_trimmed_1.fastq.gz \
        -2 02_TRIIMMED_DATA//SRR11674041_trimmed_2.fastq.gz \
        -S 04_POLISH/assembly.sam \
        -x 04_POLISH//assembly \
        --threads $SLURM_CPUS_PER_TASK \
        --no-unal

module purge 
module load samtools/1.6.1

samtools view -F 4 -bS 04_POLISH/assembly.sam | samtools sort > 04_POLISH/assembly.bam
samtools index 04_POLISH/assembly.bam

module purge
module load biojava/1.8

java -Xmx128G -jar /scratch/project_2001499/envs/pilon/pilon-1.24.jar \
        --genome 03_ASSEMBLY/assembly.fasta \
        --bam 04_POLISH/assembly.bam \
        --outdir 04_POLISH/ \
	--output pilon \
        --threads $SLURM_CPUS_PER_TASK \
        --changes
        --output pilon \
        --threads $SLURM_CPUS_PER_TASK \
        --changes