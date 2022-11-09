# Practicals

__Table of Contents:__
1. [Setting up](#setting-up-the-course-folders)
2. [Interactive use of Puhti](#interactive-use-of-puhti)
3. [QC and trimming for Nanopore reads](#qc-and-trimming-for-nanopore-reads)
4. [Metagenomic assembly with metaFlye](#metagenomic-assembly-with-metaflye)
4. [QC and trimming for Illumina reads](#qc-and-trimming-for-illumina-reads)
5. [Genome assembly with Spades](#genome-assembly-with-spades)
6. [Eliminate contaminant contigs with Kaiju](#eliminate-contaminant-contigs-with-kaiju)
7. [Assembly QC](#assembly-qc)
8. [Calculate the genome coverage](#calculate-the-genome-coverage)
9. [Genome completeness and contamination](#genome-completeness-and-contamination)
10. [Genome annotation with Prokka](#genome-annotation-with-prokka)
11. [Name the strain](#name-the-strain)
12. [Pangenomics](#pangenomics-with-anvio)
13. [Detection  of secondary  metabolites biosynthesis gene clusters](#detection-of-secondary-metabolites-biosynthesis-gene-clusters)
14. [Comparison of secondary metabolites biosynthesis gene clusters](#comparison-of-secondary-metabolites-biosynthesis-gene-clusters)

## Setting up the course folders
The main course directory is located in `/scratch/project_2001499`.  
There you will set up your own directory where you will perform all the tasks for this course.  

First list all projects you're affiliated with in CSC.

```
csc-workspaces
```

You should see the course project `MBDP_Metagenomics_2022`.
So let's create a folder for you inside the scratch folder, you can find the path in the output from the previous command.

```bash
cd PATH/TO/COURSE/SCRATCH
mkdir $USER
```

Check with `ls`; which folder did `mkdir $USER` create?

This directory (`/scratch/project_2001499/your-user-name`) is your own working directory.  
Every time you log into Puhti, you should use `cd` to navigate to this directory.

Go to your own folder and clone the course repository there. 

```bash
git clone https://github.com/MBDP-bioinformatics-courses/MBDP_Metagenomics_2022.git
```

Check what was downloaded and go to that folder. Then again check what is inside. 
**All the scripts are to be run in this folder**.  

The raw data used on this course can be found in `/scratch/project_2001499/Data/`.  
Instead of copying the data we will use links to this folder in all of the needed tasks.  
Why don't we want 14 students copying data to their own folders?


We have both Nanopore long-read data and Illumina short-read data from one of the sewage treatment plants in the article. Make separate folders for both data types to your own data folder, `01_DATA`. 

```bash
cd 01_DATA
mkdir Nanopore
mkdir Illumina
```

Then go first to the `Nanopore` folder and make a softlink to the raw long-read data.  
After making the softlink, check the content of the folder with `ls`. 

```bash
cd Nanopore
ln -s /scratch/project_2001499/Data/Nanopore/SRR11673980.fastq.gz .
cd ..
```

Then the same for the short-read data.  
Again, after making the softlink, check the content of the folder with `ls`. 

```bash
cd Illumina 
ln -s /scratch/project_2001499/Data/Illumina/*.fastq.gz .
cd ..
```

How many files (links) did you have in the Nanopore folder? What about the Illumina?


## Interactive use of Puhti

Puhti uses a scheduling system called SLURM. Most jobs are sent to the queue, but smaller jobs can be run interactively.

Interactive session is launched with `sinteractive`   .   
You can specify the resources you need for you interactive work interactively with `sinteractive -i`. Or you can give them as options to `sinteractive`.  
You always need to specify the accounting project (`-A`, `--account`). Otherwise for small jobs you can use the default resources (see below).

| Option | Function | Default | Max |  
| --     | --       | --      | --  |  
| -i, --interactive | Set resources interactively |  |  |  
| -t,  --time | Reservation in minutes or in format d-hh:mm:ss | 24:00:00 | 7-00:00:00 |
| -m, --mem | Memory in Mb       | 2000     | 76000  |  
| -j, --jobname |Job name       | interactive     |   |  
| -c, --cores     | Number of cores       | 1      | 8  |  
| -A, --account     | Accounting project       |       |  |  
| -d, --tmp     | $TMPDIR size (in GiB)      |  32     | 760  |  
| -g, --gpu     | Number of GPUs       | 0     | 0 |  


[__Read more about interactive use of Puhti.__](https://docs.csc.fi/computing/running/interactive-usage/#sinteractive-in-puhti)   

## QC and trimming for Nanopore reads

The QC for the Nanopore reads can be done with NanoPlot. It is a plotting tools for long read sequencing data and alignments. You can read more about it from here: [NanoPlot](https://github.com/wdecoster/NanoPlot).

This run will require computing resources, you can use `sinteractive` to log in to a computing node:

```bash
sinteractive -A project_2001499 -m 50G -c 4 
```

NanoPlot is not pre-installed to Puhti, but has been installed for the course under `/scratch/project_2001499/envs/nanoQC` and can be run from there.

Generate graphs for visualization of reads quality and length distribution

```bash
/scratch/project_2001499/envs/nanoQC/bin/NanoPlot -o 01_DATA/nanoplot_out -t 4 -f png --fastq softlink-to/raw_nanopore_reads.fastq.gz
```

Transfer to the output from NanoPlot (`NanoPlot-report.html`) to your own computer and open it with any browser.

Based on the results from NanoPlot:
* How many reads do we have?
* How much data we have overall?
* How the read length distribution looks like and what is the mean read length?
* How about the read quality? 
* Can you see any pattern between read length and read quality?

### Trimming and quality filtering of reads

As in the article, we will use [Filtlong](https://github.com/rrwick/Filtlong) and [Porechop](https://github.com/rrwick/Porechop) for quality and adapter trimmming of the nanopore reads.

As we have quite a lot of data, trimming also takes a while. Submitting it as batch job means we can use more resources and the job shoud also run faster.  
Make a batch job script using the example and store the batch job script under `scripts`.  
But remember to check all file paths before you submit the job. 

```bash
#!/bin/bash -l
#SBATCH --job-name trimming
#SBATCH --output 00_LOGS/trimming_out_%j.txt
#SBATCH --error 00_LOGS/trimming_err_%j.txt
#SBATCH --time 24:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 50G
#SBATCH --account project_2001499
#SBATCH --gres=nvme:200

/scratch/project_2001499/envs/nanoFiltering/bin/filtlong \
        --min_length 4000 \
        --min_mean_q 80 soflink-to/your_raw_nanopore_reads.fastq.gz |\
         gzip > 02_TRIMMED_DATA/SRR11673980_trimmed.fastq.gz

/scratch/project_2001499/envs/nanoFiltering/bin/porechop \
        -i 02_TRIMMED_DATA/SRR11673980_trimmed.fastq.gz \
        -o 02_TRIMMED_DATA/SRR11673980_chop.fastq.gz \
        --threads $SLURM_CPUS_PER_TASK
```

Wheen you are sure everything is ok, submit the job with `sbatch`. 

```bash
sbatch scripts/YOUR_SCRIPT.sh
```

After you have submitted the job, the system will print you the job id. Write that down.  
You can also check the status of your job with `squeue`. 

```bash
squeue -l -u $USER
```

And when the job is finished, you can use `seff` to see how much resources were actually used.

```bash
seff JOBID
```

After the trimming has finished and everythin looks ok, we can move on. 

### Visualizing the trimmed data with NanoPlot

It is always good idea to check that the trimming step did what it was supposed to do. So we'll the QC step on the trimemd data.   

```bash
NanoPlot -o 02_TRIMMED_DATA/NANOPLOT -t 4 -f png --fastq 02_TRIMMED_DATA/SRR11673980_chop.fastq.gz
```

And if it looks good as well, we can move on to the assembly step.  

## Metagenomic assembly with metaFlye

Fo the assembly we will use [Flye](https://github.com/fenderglass/Flye). Flye is a long-read de novo assembler that handles also metagenomic data.  
Batch job script for assembly with metaFlye:

```bash
#!/bin/bash -l
#SBATCH --job-name assembly
#SBATCH --output 00_LOGS/assembly_out_%j.txt
#SBATCH --error 00_LOGS/assembly_err_%j.txt
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

```

## QC and trimming for Illumina reads

QC for the raw data takes few minutes, depending on the allocation.  
QC does not require lot of memory and can be run on the interactive nodes using `sinteractive`.

Open interactive node and load the biokit module:

```bash
sinteractive -A project_2001499 -c 4
module load biokit
```

### Running FastQC and MultiQC

Run FastQC with the files stored in the Illumina folder. What does the `-o` and `-t` flags refer to?

```bash
mkdir 01_DATA/FASTQC
fastqc -o 01_DATA/FASTQC 01_DATA/Illumina/*.fastq.gz -t 4
```

FastQC makes a report for each file and as we have several, we can use MultiQC to combine them into one report.   
MultiQC is not included in biokit, so we need to load it separately. It is wise to first unload all other modules with `module purge`. 

```bash
module purge
module load multiqc/1.12  

multiqc --interactive -o 01_DATA/ 01_DATA/FASTQC/*
```

Copy the resulting HTML file to your local machine.
Have a look at the QC report with your favourite browser.  

After inspecting the output, it should be clear that we need to do some trimming.  
__What kind of trimming do you think should be done?__

### Running Cutadapt

We will do the trimming with [cutadapt]((http://cutadapt.readthedocs.io).
Cutadapt module has to be loaded separately.

```bash
module purge
module load cutadapt/3.5 
```

Have a look at the help pages for cutadapt. 

```
cutadapt -h
```

The adapter sequences that you want to trim are located after `-a` and `-A`.  
What is the difference with `-a` and `-A`?  
And what is specified with option `-p` or `-o`?
And how about `-m` and `-j`?  
You can also find the answers from Cutadapt [manual](http://cutadapt.readthedocs.io).

Cutadapt can handle paired-end reads, but since we have several files, we will make a for loop that trims each sample separately. The for loop will replace the `${file}` in the script with each of the sample IDs specified on the first line.  

```bash
for file in SRR11674041 SRR11674042 SRR11674043
do
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o 02_TRIMMED_DATA/${file}_trimmed_1.fastq.gz \
        -p 02_TRIMMED_DATA/${file}_trimmed_2.fastq.gz \
        01_DATA/Illumina/${file}_1.fastq.gz \
        01_DATA/Illumina/${file}_2.fastq.gz \
        --minimum-length 50 \
        --cores 4 \
        > 00_LOGS/${file}_cutadapt.log
done
```

You could now check each of the cutadapt log files and answer:

* How many read pairs we had originally?
* How many reads contained adapters?
* How many read pairs were removed because they were too short?
* How many base calls were quality-trimmed?
* Overall, what is the percentage of base pairs that were kept?

### Running FastQC and MultiQC on the trimmed reads

Then make a new folder (`FASTQC`) for the QC files of the trimmed data and run fastQC and multiQC again as you did before trimming:

FastQC:
```bash
module purge
module load biokit

mkdir 02_TRIMMED_DATA/FASTQC
fastqc -o 02_TRIMMED_DATA/FASTQC 02_TRIMMED_DATA/*_trimmed_?.fastq.gz  -t 4
```

MultiQC:
```bash
module purge
module load multiqc/1.12  

multiqc --interactive -o 02_TRIMMED_DATA/ 02_TRIMMED_DATA/FASTQC/*
```

Copy the resulting HTML file to your local machine as earlier and look how well the trimming went.  
Did you find problems with the sequences? We can further proceed to quality control using Prinseq.

To leave the interactive node, type `exit`.  


## Polishig long-read assembly with short-read data

Mapping and pilon steps here. __Still draft.__

```bash
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
``` 

## Assembly QC

After the metagenomic assembly has finished we will use the metagenomic version of Quality Assessment Tool for Genome Assemblies, [Quast](http://quast.sourceforge.net/) for (comparing and) evaluating our assembly. 

```bash
module purge
module load quast/5.2.0 
```

```bash
metaquast.py \
    -o 03_ASSEMBLY/QUAST \
    03_ASSEMBLY/assembly.fasta \
    --fast \
    --threads 4
```

Now you can move the file ` 03_ASSEMBLY/QUAST/report.html` to your computer and look for the quality control files in the web browser of your preference.


__FOR LATER USE:__
## Genome completeness and contamination

Now we have calculated different metrics for our genomes, but they don't really tell anything about the "real" quality of our genome.  
We will use checkM to calculate the completeness and possible contamination in our genomes.  
Allocate some resources (>40G memory & 4 threads) and run checkM (v. 1.1.3.) from a singularity container.  

Before running checkM, it might be good to put all genomes to one folder.

```
singularity exec --bind $PWD:$PWD,$TMPDIR:/tmp /projappl/project_2005590/containers/checkM_1.1.3.sif \
              checkm lineage_wf -x fasta PATH/TO/GENOME/FOLDER OUTPUT/FOLDER -t 4 --tmpdir /tmp
```

If you missed the output of checkM, you can re-run just the last part with:

```
singularity exec --bind $PWD:$PWD,$TMPDIR:/tmp /projappl/project_2005590/containers/checkM_1.1.3.sif \
              checkm qa ./OUTPUT/lineage.ms ./OUTPUT
```

## Genome annotation with Bakta 

__ANTTI WILL ADD PARTS HERE AND COPY THE SOFTWARE AND DB TO THE COURSE FOLDER__

Now we can annotate some of our MAGs using [Bakta](https://github.com/oschwengers/bakta).  

Allocate some resources for the job. 

```bash
sinteractive -A project_2006616 -c 4
```

And then run Bakta on your favourite MAGs.  

```bash
/scratch/project_2001499/envs/bakta/bin/bakta \
       INPUT \
       --db /scratch/project_2001499/databases/bakta/db/ \
       --prefix GENOME_NAME \
       --threads 4 \
       --output OUTPUT
```

## Metaphlan JENNI IS WORKING ON THIS

Next we will also analyze individual reads in addition to the assembly based approaches. Which files would you use for this? 

We will use a tool called [Metaphlan4](https://github.com/biobakery/biobakery/wiki/metaphlan4) to analyze these reads. 

# run gtdbtk
```bash 
export GTDBTK_DATA_PATH=/scratch/project_2005590/databases/GTDB/release202/

singularity exec --bind $GTDBTK_DATA_PATH:$GTDBTK_DATA_PATH,$PWD:$PWD,$TMPDIR:/tmp  \
                    /projappl/project_2005590/containers/gtdbtk_1.7.0.sif \
                    gtdbtk classify_wf -x fasta --genome_dir PATH/TO/GENOME/FOLDER \
                    --out_dir OUTPUT/FOLDER --cpus 4 --tmpdir gtdb_test

```

## Anvi'o

```bash 
sinteractive -A project_2006616 --cores 6 --mem 50G --tmp 100
module load anvio/7.1
```

```bash
anvi-script-reformat-fasta \
    --min-len 1000 \
    --simplify-names \
    -o 05_ANVIO/contigs.fasta \
    --report-file 05_ANVIO/reformat_report.txt \
    04_POLISH/INPUT.fasta

anvi-gen-contigs-database \
    -f 05_ANVIO/contigs.fasta \
    -o 05_ANVIO/CONTIGS.db \
    -n Long-read assembly \
    -T 6

anvi-run-hmms \
    -c 05_ANVIO/CONTIGS.db \
    -T 6

anvi-run-ncbi-cogs \
    -c 05_ANVIO/CONTIGS.db \
    -T 6

anvi-run-scg-taxonomy \
    -c 05_ANVIO/CONTIGS.db \
    -T 6

anvi-estimate-scg-taxonomy \
    -c 05_ANVIO/CONTIGS.db \
    --metagenome-mode
```

Mapping script:

```bash
#!/bin/bash
module load bowtie2/2.4.4  
bowtie2-build 05_ANVIO/contigs 05_ANVIO/contigs.fasta

for file in SRR11674041 SRR11674042 SRR11674043
do
    module load bowtie2/2.4.4
    bowtie2 \
        -1 02_TRIMMED_DATA/${file}_trimmed_1.fastq.gz \
        -2 02_TRIMMED_DATA/${file}_trimmed_2.fastq.gz \
        -x 05_ANVIO/contigs
        -S 05_ANVIO/${file}.sam \
        --threads 6 \
        --no-unal

    module purge
    module load samtools

    samtools view -F 4 -bS 05_ANVIO/${file}.sam |\
        samtools sort > 05_ANVIO/${file}.bam
    
    samtools index 05_ANVIO/${file}.bam
    
    anvi-profile-blitz \
        -i 05_ANVIO/${file}.bam \
        -c 05_ANVIO/CONTIGS.db \
        -S ${file} \
        -o 05_ANVIO/${file}_PROFILE \
        -T 6
done 
```

Run script

```bash
bash scripts/map2assembly.sh
```

```bash
anvi-merge \
    -o 05_ANVIO/SAMPLES-MERGED \
    -c 05_ANVIO/CONTIGS.db \
    05_ANVIO/*_PROFILE/PROFILE.db 
```

Mini manual for screen:

    screen -S NAME - open a screen and give it a session name NAME
    screen - open new screen without specifying any name
    screen -ls - list all open sessions
    ctrl + a + d - to detach from a session (from inside the screen)
    screen -r NAME - re-attach to a detached session using the name
    screen -rD - re-attach to a attached session
    exit - close the screen and kill all processes running inside the screen (from inside the screen)


Then you can log out and log in again, but this time in a bit different way.
You need to specify your PORT and the NODEID to which you connected and also the NUMBER of the login node you where your screen is running. Also change your username in the command below.

```
ssh -L PORT:NODEID.bullx:PORT USERNAME@puhti-loginNUMBER.csc.fi
```
And in windows using Putty:
In SSH tab select "tunnels". Add:

Source port: PORT
Destination: NODEID.bullx:PORT
Click add and connect to the right login node, login1 or login2.

Then go back to your screen and launch the interactive interface.
Remember to change the PORT.
