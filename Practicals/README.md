# Practicals

__Table of Contents:__
1. [Setting up](#setting-up-the-course-folders)
2. [Connnecting to Puhti](#connecting-to-puhti-with-visual-studio-code)
3. [Interactive use of Puhti](#interactive-use-of-puhti)
4. [QC and trimming for Nanopore reads](#qc-and-trimming-for-nanopore-reads)
5. [Metagenomic assembly with metaFlye](#metagenomic-assembly-with-metaflye)
6. [QC and trimming for Illumina reads](#qc-and-trimming-for-illumina-reads)
7. [Assembly polishing](#polishing-long-read-assembly-with-short-read-data)
8. [Assembly QC](#assembly-qc)
9. [Taxonomic composition](#taxonomic-composition-with-metaphlan-using-short-reads)
10. [Genome-resolved metagenomics](#genome-resolved-metagenomics-with-anvio)
11. [Genome taxonomy](#taxonomic-annotation-of-mags-with-gtdb-tk)
12. [Genome annotation with Bakta](#genome-annotation-of-mags-with-bakta)
13. [Viromics](#identifying-viral-contings-from-the-metagenome)
14. [Visualizing Metaphlan results](#visualizing-metaphlan-results-with-r)

## Connecting to Puhti with Visual Studio Code

* Launch Visual Studio Code
* _Only on the first time:_ 
    - _Open `Extensions` (one of the icons on the left) and install `Remote - SSH`_  
* Down left corner you will have a (green) button with "><" (hoover over it and it says "Open a Remote Window"), click it 
* Choose "Connect Current Window to Host..."
* Type in the **user<span>@puhti.csc.fi** and hit "Enter" (change "user" for your own CSC username) 
* Type your password and hit "Enter"
* In the following dialogue, type **yes** and hit "Enter"

When the down left corner says `SSH:puhti.csc.fi`, you're connected.
* From the menu select `Terminal > New Terminal` and you should see a new panel. This is the __command line__.

* When you need to logout just type **exit** in the terminal/command line and hit "Enter"  
(or you can click the down left corner and choose "Close Remote Connection")

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

Mini manual for screen: 
- `screen -S NAME` - open a screen and give it a session name NAME
- `screen` - open new screen without specifying any name
- `screen -ls` - list all open sessions
- `ctrl + a + d` - to detach from a session (from inside the screen)
- `screen -r NAME` - re-attach to a detached session using the name
- `screen -rD` - re-attach to a attached session
- `exit` - close the screen and kill all processes running inside the screen (from inside the screen)

## QC and trimming for Nanopore reads

The QC for the Nanopore reads can be done with NanoPlot. It is a plotting tools for long read sequencing data and alignments. You can read more about it from here: [NanoPlot](https://github.com/wdecoster/NanoPlot).

This run will require computing resources, you can use `sinteractive` to log in to a computing node:

```bash
sinteractive -A project_2001499 -m 50G -c 4 
```

NanoPlot is not pre-installed to Puhti, but has been installed for the course under `/scratch/project_2001499/envs/nanoQC` and can be run from there.

Generate graphs for visualization of reads quality and length distribution

```bash
/scratch/project_2001499/envs/nanoQC/bin/NanoPlot \
    -o 01_DATA/nanoplot_out \
    -t 4 \
    -f png \
    --fastq softlink-to/raw_nanopore_reads.fastq.gz
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
#SBATCH --output 00_LOGS/trimming-%j.out
#SBATCH --error 00_LOGS/trimming-%j.err
#SBATCH --time 24:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 80G
#SBATCH --account project_2001499
#SBATCH --gres=nvme:200

/scratch/project_2001499/envs/nanoFiltering/bin/filtlong \
        --min_length 4000 \
        --min_mean_q 80 \
        soflink-to/your_raw_nanopore_reads.fastq.gz |\
         gzip > 02_TRIMMED_DATA/SRR11673980_trimmed.fastq.gz

/scratch/project_2001499/envs/nanoFiltering/bin/porechop \
        -i 02_TRIMMED_DATA/SRR11673980_trimmed.fastq.gz \
        -o 02_TRIMMED_DATA/SRR11673980_chop.fastq.gz \
        --threads $SLURM_CPUS_PER_TASK
```

When you are sure everything is ok, submit the job with `sbatch`. 

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
sinteractive -A project_2001499 -m 50G -c 4 
```

```bash
/scratch/project_2001499/envs/nanoQC/bin/NanoPlot \
    -o 02_TRIMMED_DATA/NANOPLOT \
    -t 4 \
    -f png \
    --fastq 02_TRIMMED_DATA/SRR11673980_chop.fastq.gz
```

And if it looks good as well, we can move on to the assembly step.  

## Metagenomic assembly with metaFlye

Fo the assembly we will use [Flye](https://github.com/fenderglass/Flye). Flye is a long-read de novo assembler that handles also metagenomic data.  
Before you start the assembly, have a look at the Flye manual, escpecially the parts about Nanopore data and metagenome assembly.  

After you have read the manual you can make a batch job script for Flye in the `scripts` folder. 

Batch job script for assembly with metaFlye:

```bash
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
```

Submit the job to the queue.

```bash 
sbatch scripts/YOUR_SCRIPT_NAME
```

Using the commands from above, you can check if your job is already running.  

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

We will do the trimming with [cutadapt](http://cutadapt.readthedocs.io).
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
        -o 02_TRIMMED_DATA/${file}_trimmed_R1.fastq.gz \
        -p 02_TRIMMED_DATA/${file}_trimmed_R2.fastq.gz \
        01_DATA/Illumina/${file}_1.fastq.gz \
        01_DATA/Illumina/${file}_2.fastq.gz \
        --minimum-length 50 \
        --cores 4 \
        > 00_LOGS/cutadapt_${file}.log
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
fastqc -o 02_TRIMMED_DATA/FASTQC 02_TRIMMED_DATA/*_trimmed_R?.fastq.gz  -t 4
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


## Polishing long-read assembly with short-read data

It is possible to polish the long-read assembly using either the original reads, as was done in the article, or by using better quality data, such as short-reads, from the same sample.  

However, the polishing step is very slow, so we won't run it on the course. We will rely on the polishing step done by metaFlye. 

But as an example, here are the steps for polishing our assembly using the trimmed short-read data from the same sample.

```bash
# module load bowtie2/2.4.4 
# bowtie2-build 03_ASSEMBLY/assembly.fasta \
#                 04_POLISH/assembly

# bowtie2 -1 02_TRIMMED_DATA/SRR11674041_trimmed_R1.fastq.gz \
#         -2 02_TRIMMED_DATA/SRR11674041_trimmed_R2.fastq.gz \
#         -S 04_POLISH/assembly.sam \
#         -x 04_POLISH/assembly \
#         --threads $SLURM_CPUS_PER_TASK \
#         --no-unal

# module purge 
# module load samtools/1.6.1

# samtools view -F 4 -bS 04_POLISH/assembly.sam | samtools sort > 04_POLISH/assembly.bam
# samtools index 04_POLISH/assembly.bam

# module purge
# module load biojava/1.8

# java -Xmx128G -jar /scratch/project_2001499/envs/pilon/pilon-1.24.jar \
#         --genome 03_ASSEMBLY/assembly.fasta \
#         --bam 04_POLISH/assembly.bam \
#         --outdir 04_POLISH/ \
#         --output pilon \
#         --changes
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
    --max-ref-number 0 \
    --fast \
    --threads 4
```

Now you can move the file ` 03_ASSEMBLY/QUAST/report.html` to your computer and look for the quality control files in the web browser of your preference.  

## Taxonomic composition with Metaphlan using short-reads 

Next we will also analyze individual reads in addition to the assembly based approaches. Which files would you use for this? 

We will use a tool called [Metaphlan4](https://github.com/biobakery/biobakery/wiki/metaphlan4) to analyze these reads. From the paired end Illumina reads only R1 reads are used for the following analyses. Metaphlan takes fastq or fasta-formatted files but cannot read compressed files. Our filtered data is in fastq.gz format so it we need to do some uncompressing with command `gunzip`. 

Use flag `-h` for help (or Google) to check how the command works and uncompress your files.

```bash
#!/bin/bash -l
#SBATCH --job-name metaphlan
#SBATCH --output 00_LOGS/metaphlan-%j.out
#SBATCH --error 00_LOGS/metaphlan-%j.err
#SBATCH --time 06:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 100G
#SBATCH --account project_2001499


# load metaphlan4

module load metaphlan/4.0.2

# compare the three fastq-files against metaphlan databases we have downloaded and formatted to the databases-folder

for file in SRR11674041 SRR11674042 SRR11674043
do
    metaphlan \
       02_TRIMMED_DATA/${file}_trimmed_R1.fastq \
        --input_type fastq \
        --bowtie2db /scratch/project_2001499/databases/metaphlan/ \
        --bowtie2out 07_METAPHLAN/${file}.bowtie2.bz2 \
        --nproc $SLURM_CPUS_PER_TASK \
        -o 07_METAPHLAN/${file}_metaphlan.txt \
      
done
```    

## Genome-resolved metagenomics with anvi'o

Anvi'o is an analysis and visualization platform for omics data. We will use anvi'o for binning contigs into metagenome-assembled genomes (MAGs).  
You should definitely take a look at their [website](https://anvio.org/) and maybe even join their [discord channel](https://discord.gg/C6He6mSNY4).


### Reformat the assembly file 

Before creating the contigs database in anvi'o, we need to do some reformatting for our assmebly file.  
The program removes contigs shorter than 1 000 bp and simplifyes the sequence names. 

```bash 
module load anvio/7.1

anvi-script-reformat-fasta \
    --min-len 1000 \
    --simplify-names \
    -o 05_ANVIO/contigs.fasta \
    --report-file 05_ANVIO/reformat_report.txt \
    03_ASSEMBLY/assembly.fasta
```

### Nucleotide composition and the contigs database

The first step in our genome-resolved metagenomics pipeline is the contruction of contigs database from our metagenomic contigs. During this step anvi'o calls genes, calculates the nucleotide composition for each contigs and annotates the identified genes in three different steps.  

We will run the first steps as a batch job. 

Batch job script. Save it to the `scripts` folder. 

```bash
#!/bin/bash -l
#SBATCH --job-name anvi-contigs-db
#SBATCH --output 00_LOGS/anvi-contigs-db-%j.out
#SBATCH --error 00_LOGS/anvi-contigs-db-%j.err
#SBATCH --time 03:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 30G
#SBATCH --account project_2001499
#SBATCH --gres=nvme:200

module purge
module load anvio/7.1

# generate contigs DB
anvi-gen-contigs-database \
    -f 05_ANVIO/contigs.fasta \
    -o 05_ANVIO/CONTIGS.db \
    -n LongReadAssembly \
    -T $SLURM_CPUS_PER_TASK

# annotate marker genes
anvi-run-hmms \
    -c 05_ANVIO/CONTIGS.db \
    -T $SLURM_CPUS_PER_TASK

# annotate COGs
anvi-run-ncbi-cogs \
    -c 05_ANVIO/CONTIGS.db \
    --cog-data-dir /scratch/project_2001499/databases/anvio/ \
    -T $SLURM_CPUS_PER_TASK

# annotate single-copy core genes
anvi-run-scg-taxonomy \
    -c 05_ANVIO/CONTIGS.db \
    --scgs-taxonomy-data-dir /scratch/project_2001499/databases/anvio/ \
    -T $SLURM_CPUS_PER_TASK
```

Submit the job.

``` bash 
sbatch scripts/YOUR_SCRIPT_NAME
```

### Differential coverage by mapping

The differential coverage for each contig is calculated by mapping sequencing reads to the assembly.  
We will use Bowtie2 to map the short-read Illumina data to our assembly (the anvio-reformatted version of it). 

Mapping batch job script, again save it to the `scripts` folder:

```bash
#!/bin/bash -l
#SBATCH --job-name anvi-mapping
#SBATCH --output 00_LOGS/anvi-mapping-%j.out
#SBATCH --error 00_LOGS/anvi-mapping-%j.err
#SBATCH --time 2:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 40
#SBATCH --mem 80G
#SBATCH --account project_2001499
#SBATCH --gres=nvme:200

module load bowtie2/2.4.4  
bowtie2-build --threads $SLURM_CPUS_PER_TASK 05_ANVIO/contigs.fasta 05_ANVIO/contigs

for file in SRR11674041 SRR11674042 SRR11674043
do
    module load bowtie2/2.4.4
   
   bowtie2 \
        -1 02_TRIMMED_DATA/${file}_trimmed_R1.fastq.gz \
        -2 02_TRIMMED_DATA/${file}_trimmed_R2.fastq.gz \
        -x 05_ANVIO/contigs \
        -S 05_ANVIO/${file}.sam \
        --threads $SLURM_CPUS_PER_TASK \
        --no-unal

    module purge
    module load samtools

    samtools view -@ $SLURM_CPUS_PER_TASK -F 4 -bS 05_ANVIO/${file}.sam |\
        samtools sort -@ $SLURM_CPUS_PER_TASK > 05_ANVIO/${file}.bam
    
    samtools index -@ $SLURM_CPUS_PER_TASK 05_ANVIO/${file}.bam
    
    rm 05_ANVIO/${file}.sam
done 
```

Submit the job.

```bash
sbatch scripts/YOUR_SCRIPT_NAME
```

### Profile databases and merging the profile databases

When the contigs database and all three mappings are ready, we can make the profile databases and finally merge all three profile databases into one. 

Once again, save the batch job script to the `scripts` folder.

```bash
#!/bin/bash -l
#SBATCH --job-name anvi-profiling
#SBATCH --output 00_LOGS/anvi-profiling-%j.out
#SBATCH --error 00_LOGS/anvi-profiling-%j.err
#SBATCH --time 2:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 80G
#SBATCH --account project_2001499
#SBATCH --gres=nvme:200

for file in SRR11674041 SRR11674042 SRR11674043
do
    module load anvio/7.1
    
    anvi-profile \
        -i 05_ANVIO/${file}.bam \
        -c 05_ANVIO/CONTIGS.db \
        -S ${file} \
        --min-contig-length 5000 \
        -o 05_ANVIO/${file}_PROFILE \
        -T $SLURM_CPUS_PER_TASK
done

anvi-merge \
    -o 05_ANVIO/SAMPLES-MERGED \
    -c 05_ANVIO/CONTIGS.db \
    --enforce-hierarchical-clustering \
    05_ANVIO/*_PROFILE/PROFILE.db 
```

Submit the job. 

```bash
bash scripts/YOUR_SCRIPT_NAME
```

### Quick look at the taxonomic content of our metagenome

While the profile job is running, we can quickly estimate the taxonomic composition of our assembled metagenome based on the single-copy cores genes found in there. 

```bash
module load anvio/7.1

anvi-estimate-scg-taxonomy \
    -c 05_ANVIO/CONTIGS.db \
    --metagenome-mode  \
    --output-file 05_ANVIO/scg-taxonomy.txt
```

Have a look at the output.

```bash 
less -S 05_ANVIO/scg-taxonomy.txt
```

### Manual binning and the anvi'o interactive interface

Then everything should  be ready and we can start the manual binning part.  
We will go thru the first steps together. And provide you the port number that you will need in the next steps. 

```bash 
sinteractive -A project_2001499 --mem 20G
module load anvio/7.1
export $ANVIOPORT=YOUR_PORT_HERE
```

Start interactive interface

```bash
anvi-interactive \
    -c 05_ANVIO/CONTIGS.db \
    -p 05_ANVIO/SAMPLES_MERGED/PROFILE.db \
    --server-only \
    --port-number $ANVIOPORT
```

After making the pre-clustering, start refining the clusters.

```bash
anvi-refine \
    -c 05_ANVIO/CONTIGS.db \
    -p 05_ANVIO/SAMPLES-MERGED/PROFILE.db \
    --server-only \
    --port-number $ANVIOPORT \
    --collection-name PreCluster \
    --bin-id Bin_1
```

When all cluster have been refined a bit more, we can make a new collection from these. 

```bash
anvi-rename-bins \
    -c 05_ANVIO/CONTIGS.db \
    -p 05_ANVIO/SAMPLES-MERGED/PROFILE.db \
    --collection-to-read PreCluster \
    --collection-to-write PreBins \
    --prefix Preliminary \
    --report-file 05_ANVIO/PreBin_report.txt
```
And then we can make a summary of the cluster or bins we have so far. 

```bash
anvi-summarize \
    -c 05_ANVIO/CONTIGS.db \
    -p 05_ANVIO/SAMPLES-MERGED/PROFILE.db \
    --collection-name PreBins \
    --output-dir 05_ANVIO/PreBins_SUMMARY \
    --quick-summary
```

Download the output folder to your own computer and have a look at the `index.html` file in the output folder.  
__Do you have any promising bins?__  
Then proceed with refining the most promising bins that you have with `anvi-refine`. Remember to specify the new collection name and check the new bin names.  

The goal is to get at least few good bins that we can call MAGs. 

Then when you're happy with your work, you can run the `anvi-summarize` command on your collection without the `--quick-summary` option. 

### Useful commands for binning workflow

__Refine a bin__

```bash
anvi-refine
```

__Rename all bins and make a new collection from those__

```bash
anvi-rename-bins
```

__Summarize a collection__

```bash
anvi-summarize 
```

## Taxonomic annotation of MAGs with GTDB-Tk

After binning you hopefully have some good quality MAGs. The next step is to give some names to these metagenome assembled genomes. We will use Genome Taxonomy Database ([GTDB](https://gtdb.ecogenomic.org/)) and a tool called [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/installing/index.html#installing-gtdbtk-reference-data) for this.   

GTDB-Tk needs a database for the taxonomic annotation and that has been already done. 
Instructions are found from [here](https://ecogenomics.github.io/GTDBTk/installing/index.html#installing-gtdbtk-reference-data).  

And the actual commands were:

```bash
#############  (THIS HAS BEEN DONE ALREADY)  #########################
## download gtdb database
# wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz
# tar -xzf gtdbtk_v2_data.tar.gz
#############  (THIS HAS BEEN DONE ALREADY)  #########################
```

If you're connected to a interactive node, exit from that with `exit`. We need more resources for the next step. 

```bash
sinteractive -A project_2001499 --mem 70G -c 6
```
First we need to store to location of the database to `$GTDBTK_DATA_PATH` environment variable. Otherwise GTDB-Tk wonät find it and will complain. 

```bash
export GTDBTK_DATA_PATH=/scratch/project_2001499/databases/GTDB/releaseXX
```

Then we can have a look at the different options of the tool.

```bash 
singularity exec --bind --bind $GTDBTK_DATA_PATH:$GTDBTK_DATA_PATH,$PWD:$PWD,$TMPDIR:/tmp \
    /scratch/project_2001499/envs/gtdb-tk.sif \
    gtdbtk -h
```

And then actually run the program with few selected MAGs. The MAG sequences should all be in one folder. 

``` bash
singularity exec --bind --bind $GTDBTK_DATA_PATH:$GTDBTK_DATA_PATH,$PWD:$PWD,$TMPDIR:/tmp \
    /scratch/project_2001499/envs/gtdb-tk.sif \
    gtdbtk classify_wf \
    -x fasta \
    --genome_dir PATH/TO/GENOME/FOLDER \
    --out_dir OUTPUT/FOLDER \
    --cpus 6 \
    --tmpdir gtdb_test
```

## Genome annotation of MAGs with Bakta 

Now we can annotate some of our MAGs using [Bakta](https://github.com/oschwengers/bakta).  

Allocate some resources for the job. 

```bash
sinteractive -A project_2001499 -c 4
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

## Identifying viral contings from the metagenome

The assembled bulk metagenome contains also viral sequences. There are various bioinformatics tools for detecting viral sequences in metagenomes, and they are based on different algorithms and thus, perform differently. During this course, we will use [Virsorter2](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y) and [Lazypipe](https://www.helsinki.fi/en/projects/lazypipe). We will also look at the [What-the-Phage](https://www.biorxiv.org/content/10.1101/2020.07.24.219899v3.full) pipeline sample results (you will get ready output files).

Questions to think about:

* Why does the metagenome contain also viral DNA? 
* Which types of viruses can be seen from a bulk metagenome obtained from the environmental DNA sample? Can we find, e.g., RNA viruses there?

### Running Virsorter2

Make a directory called `Virsorter2` in your `06_VIROMICS` directory. Make a batch job script there.  
Sample script:

```bash
#!/bin/bash
#SBATCH --job-name=virsorter2
#SBATCH --account=project_2001499
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=15G
#SBATCH --partition=small
#SBATCH --error=virsorter2_err_%j.txt

export INPUT_PATH=/scratch/project_2001499/$USER/MBDP_Metagenomics_2022/03_ASSEMBLY/

apptainer exec \
     --bind $INPUT_PATH:$INPUT_PATH,$PWD:$PWD \
     /scratch/project_2001499/envs/virsorter2/virsorter2.sif \
     virsorter run \
     -w virsorter2.out \
     -i $INPUT_PATH/assembly.fasta \
     --min-length 1500 \
     --include-groups "dsDNAphage,ssDNA,RNA,NCLDV,lavidaviridae" \
     -j 30 \
     all
```
where, `virsorter.out` is the output (results) directory, you can name it yourself, and the `assembly.fasta` is the assembly input file. The option `--include-groups` specifies virus groups to search for (Virsorter version 2.2.3, which we have installed in Puhti, has only dsDNAphage and ssDNA virus groups to search for by default, we can add more). Explore more about the options from the [manual](https://github.com/jiarong/VirSorter2) or by calling 
```
apptainer exec --bind $PWD:$PWD scratch/project_2001499/envs/virsorter2/virsorter2.sif virsorter run -h
```
Note that you need to specify the path to your assembly file.

Submit a batch job.

### Virsorter2 output files (results)

The useful output files: `final-viral-boundary.tsv`, `final-viral-score.tsv`, and `final-viral-combined.fa`. Copy them to your computer and explore.

**final-viral-boundary.tsv** file contains information about contigs that were identified as viral. Full sequences identified as viral have suffix "||full"; partial sequences identified as viral have suffix "||{i}index_partial" ("{I}" can be numbers starting from 0 to max number of viral fragments found in that contig).
* How many contigs were annotated as viral? How many full and partial?
* Which virus groups were predicted? 
* The group column is the classifier of viral group that gives a high score, but can it be used as a reliable classification?
* Which groups can Virsorter2 identify (check the [publication](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y)) and which would you expect to see in our sample? Compare this to the ones predicted.

**final-viral-score.tsv** shows the score of each viral sequences across groups. The score is ranging from 0 to 1, and higher means more like to be viral. 
* Check the score for different contigs, how many have 1.000? 
* Check the number of hallmark genes found for different contigs: are there many viral contings with 0 viral genes and what is the maximum number found?
* Why some contigs with 0 viral genes are still annotated as viral? Check the Virsorter2 work principles from the [publication](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y).
* What is the length of the contigs? The smallest vs the biggest? What is the average bacteriophage genome length in principle?

**final-viral-combined.fa** contains the viral contigs sequences in the fasta format. They can be used for the downstream analyses. 

### Checking the quality of Virsorter2 predicted contigs with CheckV.

We will use [CheckV](https://www.nature.com/articles/s41587-020-00774-7) to assess the quality and completeness of the obtained viral contigs.

CheckV needs the database, it is in `/scratch/project_2001499/checkv-db` (and this is specified in the script below).

Make a `CheckV` directory in your `06_VIROMICS` directory and run CheckV interactively from there:

```bash
sinteractive -A project_2001499 -m 10G -c 8 
```
Run the job:

```bash
cd scratch/project_2001499/YOUR_VIRSORTER2_OUTPUT_DIRECTORY

apptainer exec --bind $PWD:$PWD,/scratch/project_2001499/checkv-db/db:/db \
    /scratch/project_2001499/envs/checkV/checkV.sif checkv end_to_end PATH_TO_VIRSORTER2_OUTPUT/final-viral-combined.fa \
    checkv_out -t 8 -d /db
```
Note that you need to specify the path to your Virsorter2 results directory and name the CheckV output directory (`checkv_out` in this sample script).

### CheckV output files (results)

The useful files: `quality_summary.tsv`, `completeness.tsv`, and `complete_genomes.tsv`. Copy them to your computer and explore.

**quality_summary.tsv** is the main output file.
* Are there any proviruses predicted?
* Are most contigs of low, medium or high quality? See how the quality corresponds to completeness.
* Are there 100% complete viral genomes listed?
* Are there contigs with kmer_freq > 1? This indicates that the viral genome is represented multiple times in the contig, which is quite rare.
* What warnings does the table contain?

**completeness.tsv** shows how the completeness was estimated, i.e. AAI- and HMM-based completeness. See the [publication](https://www.nature.com/articles/s41587-020-00774-7) and [FAQs](https://bitbucket.org/berkeleylab/checkv/src/master/) for more details about these two methods.

**complete_genomes.tsv** shows all viral genomes identified as complete ones.
* How many complete genomes are there in your dataset?
* What are the confidence levels for complete genomes predictions?

In practice, if you continue with downstream applications (not during this course), you might want to get rid of possible false positives in the final set of contigs identified as viral. For example, a subset of viral contigs that have at least 1 viral gene or at least 10 kbp long and 50% complete can be selected. The exact selection criteria would depend on specific research questions you have and thus, on the type of analysis you would like to perform. [Here](https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3) you can find a useful protocol for manual curation of Virsorter2 output.

### Lazypipe

[Lazypipe](https://www.helsinki.fi/en/projects/lazypipe) is a pipeline for identifying virus sequences in metagenomes, developed at the University of Helsinki. It is available as a [preinstalled module](https://docs.csc.fi/apps/lazypipe/) in Puhti. The input for the pipeline is Next Generation Sequencing data. We will use Illumina reads from the same type of samples as the long reads data that you used for the assembly with Metaflye (and Virsorter2).

Make a directory called `Lazypipe` in your `06_VIROMICS` directory and call for the Lazypipe batch job from there. Note that you need to specify the path to the Illumina reads. 

```bash
module load r-env-singularity
module load biokit
module load lazypipe

cp /appl/soft/bio/lazypipe/2.0/lazypipe/default.config.yaml config.yaml

lazypipe.pl
sbatch-lazypipe -1 /PATH_TO_ILLUMINA_READS/SRR11674042_R1.fastq.gz

```
You will be interactively asked for information that is needed to construct a batch job:
accounting project:
* maximum duration of the job (default 24 hours) - choose 24 hours
* memory reservation (default 8G) - choose 20G
* number of computing cores to use (default 8)
* email notifications - provide your email address to be notified when the batch job is finished.

### Lazypipe output files (results)

Useful output files: `qc.readsurv.jpeg`, `abund_table.xlsx`, `krona_graph.html`, and `contigs.annot.xlsx`. Download them to your computer and explore.

**qc.readsurv.jpeg** is the quality control (QC) plot tracking retained reads after each pipeline step.

**abund_table.xlsx** shows which taxa were identified in the sample and how many read pairs (readn) and contigs (consign) were assigned to some taxon. Search for viruses.
* What viral taxa can be found? 
* Do they have many contigs assigned?

**krona_graph.html** displays estimated taxon abundancies. This graph is interactive: you can zoom into different groups and see their abundance as %. Select the viruses piece of the graph to see which virus groups were identified.
* What is the abundance of viruses in the studied sample?
* Which of the identified virus groups are the most and the least abundant?
* Which bacterial groups are the most abundant?

**contigs.annot.xlsx** contains all info for the classified contigs separately for viruses (not bacteriophages), bacteria, bacteriophages, eukaryotes, and unknown. Explore the viruses and bacteriophages tabs. 
* Which putative ORF functions could be found for viral (including bacteriophages) sequences?

Viral contigs for the possible downstream analyses are found from **contigs_vi.fa**.

### What-the-Phage sample output files

“What the Phage” (WtP) is a parallel multitool pipeline for phage prediction combined with annotation and classification. It includes [11 phage-predicting tools](https://mult1fractal.github.io/wtp-documentation/#example-result-output) that are all run in parallel. The user can then compare the results of the different prediction tools summarised in charts and tables and make own decisions which dataset to proceed with. 

Due to time limits and some technical challenges associated with running this pipeline, we will only see sample output files during this course. However, you could try it later yourself, following the [authors' instructions](https://mult1fractal.github.io/wtp-documentation/installation/quick_installation/) as well as [CSC notes](https://docs.csc.fi/support/tutorials/nextflow-puhti/) about this pipeline installation. 

The sample output files for this course were obtained from the same Illumina reads dataset, as used for the Lazypipe. The reads were assembled using the [Megahit assembler](https://academic.oup.com/bioinformatics/article/31/10/1674/177884). The output summary file to download and explore: `/scratch/project_2001499/WtP_output/final_report.html`.
    
Explore different tabs of the report. It might be easier to view the tables when they are downloaded as .xlsx or .csv, rather than using the browser.

**Overview** shows a plot summarizing the prediction performance of each tool for the studied sample. Blue bars on the left show the total number of identified phage contigs per tool. Black bars show the number of contigs that each tool (or a combination of tools) has uniquely identified. A dot matrix below refers to the combinations of tools that identified the same contigs as phages.

* Which tools were applied for phage predictions?
* Which tool gave the highest number of phage contigs? Which gave the lowest? 
* Are the predicted phage contigs consistent between different tools (see the dots matrix)?

**Phage prediction table** contains a list of all predicted phage contigs and shows which tool gave which prediction and with which score.
* Do predictions from different tools overlap?

**CheckV output** lists all the predicted contigs and shows their quality based on CheckV.
* What is the typical phage contig length in this dataset? What is the lowest and highest?
* Are most contigs of high/low/medium quality?
* Are there complete contigs in the list?
* Are there proviruses?
* Are there contigs with no viral genes?

**Phage annotations** shows predicted functions obtained for phage ORFs.
* Which functions were predicted? 
* Could you try to categorise them (e.g. DNA modification, virion structural genes, etc.)?
* Which functional predictions are specific for viruses?

**Taxonomic phage classification** shows taxonomical predictions for phage contigs.
* Could any contigs be classified?

To sum up, how would you proceed with the WtP results taking into account the possibility of having false positives? For further analyses, would you
* take all predictions? Or
* take predictions only from some tools? Or
* make a subset based on the quality check? 

### Summary and downstream analyses of the identified viral sequences

Compare the data you got from __Virsorter2__, __What-the-Phage__ and __Lazypipe__. 

* What have you learnt about the viruses present in the studied datasets with each of these tools?
* What would you consider as pluses and minuses in the usage and the final data types when comparing the tested tools?
* If you need to study viruses in a metagenome in future, which of the tested tools would you select for your own work?

We are not proceeding with further viromics analyses during this course, but if you would like to test other viral identification tools and pipelines, see references in the lecture slides. 
For downstream applications, you could check, for example:
* [vConTACT v.2.0](https://pubmed.ncbi.nlm.nih.gov/31061483/) for taxonomic assignments of the identified virus genomes;
* [iVirus 2.0](https://www.nature.com/articles/s43705-021-00083-3), an analytic platform with many tools, protocols and public datasets;
* [IMG/VR](https://academic.oup.com/nar/article/49/D1/D764/5952208), the largest collection of viral sequences obtained from metagenomes.

## Visualizing Metaphlan results with R

open Rstudio
