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

You should see the course project `Metagenomics_2019`.
So let's create a folder for you inside the scratch folder, you can find the path in the output from the previous command.

```bash
cd PATH/TO/COURSE/SCRATCH
mkdir $USER
```

Check with `ls`; which folder did `mkdir $USER` create?

This directory (`/scratch/project_2001499/your-user-name`) is your working directory.  
Every time you log into Puhti, you should use `cd` to navigate to this directory, and **all the scripts are to be run in this folder**.  

__CLONE GIT REPOSITORY__ 
```bash
# code here
```

The raw data used on this course can be found in `/scratch/project_2001499/Data/`.  
Instead of copying the data we will use links to this folder in all of the needed tasks.  
Why don't we want 14 students copying data to their own folders?

__MAKE SOFTLINK__

```bash
# create directories
# make softlinks
```

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
/scratch/project_2001499/envs/nanoQC/bin/NanoPlot -o 00_DATA/nanoplot_out -t 4 -f png --fastq softlink-to/raw_nanopore_reads.fastq.gz
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

```bash
/scratch/project_2001499/envs/nanoFiltering/bin/filtlong --min_length 4000 --min_mean_q 80 soflink-to/your_raw_nanopore_reads.fastq.gz | gzip > 02_TRIMMED_DATA/SRR11673980_trimmed.fastq.gz
/scratch/project_2001499/envs/nanoFiltering/bin/porechop -i 02_TRIMMED_DATA/SRR11673980_trimmed.fastq.gz -o 02_TRIMMED_DATA/SRR11673980_chop.fastq.gz --threads 4
```

### Optional - Visualizing the trimmed data with NanoPlot
```bash
NanoPlot -o 02_TRIMMED_DATA/nanoplot_out -t 4 -f png --fastq 02_TRIMMED_DATA/SRR11673980_chop.fastq.gz
```

## Metagenomic assembly with metaFlye

__DO NOT RUN__
```bash
/scratch/project_2001499/envs/Flye/bin/flye 
        --nano-raw 02_TRIMMED_DATA/SRR11673980_chop.fastq.gz
        --threads $SLURM_CPUS_PER_TASK
        --meta
        --out-dir 03_ASSEMBLY
```

Batch job script to assemble and polish

```bash
#!/bin/bash -l
#SBATCH --job-name flye
#SBATCH --output flye_out_%j.txt
#SBATCH --error flye_err_%j.txt
#SBATCH --time 24:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 40
#SBATCH --mem 300G
#SBATCH --account project_2001499
#SBATCH --gres=nvme:200

# assembly
/scratch/project_2001499/envs/Flye/bin/flye  \
        --nano-raw 02_TRIMMED_DATA/SRR11673980_chop.fastq.gz \
        --threads $SLURM_CPUS_PER_TASK \
        --meta \
        --out-dir 03_ASSEMBLY

# polishing
/scratch/project_2001499/envs/medaka/bin/medaka_consensus \
        -i 02_TRIMMED_DATA/SRR11673980_chop.fastq.gz \
        -d 03_ASSEMBLY/assembly.fasta \
        -o medaka_out \
        -m r941_prom_high_g303 \
        --threads $SLURM_CPUS_PER_TASK
```

## QC and trimming for Illumina reads
QC for the raw data takes few minutes, depending on the allocation.  
Go to your working directory and make a folder called e.g. `fastqc_raw` for the QC reports.  

QC does not require lot of memory and can be run on the interactive nodes using `sinteractive`.

Activate the biokit environment and open interactive node:

```bash
sinteractive -A project_2001499
module load biokit
```

Now each group will work with their own sequences. Create the variables R1 and R2 to represent the path to your files. Do that just for the strain you will use:

```bash
#### Illumina Raw sequences for the cyanobacteria strain 328
R1=/scratch/project_2001499/COURSE_FILES/RAWDATA_ILLUMINA/A045-328-GGTCCATT-AGTAGGCT-Tania-Shishido-run20211223R_S45_L001_R1_001.fastq.gz
R2=/scratch/project_2005590/COURSE_FILES/RAWDATA_ILLUMINA/A045-328-GGTCCATT-AGTAGGCT-Tania-Shishido-run20211223R_S45_L001_R2_001.fastq.gz

#### Illumina Raw sequences for the cyanobacteria strain 327
R1=/scratch/project_2005590/COURSE_FILES/RAWDATA_ILLUMINA/A044-327-2-CTTGCCTC-GTTATCTC-Tania-Shishido-run20211223R_S44_L001_R1_001.fastq.gz
R2=/scratch/project_2005590/COURSE_FILES/RAWDATA_ILLUMINA/A044-327-2-CTTGCCTC-GTTATCTC-Tania-Shishido-run20211223R_S44_L001_R2_001.fastq.gz

#### Illumina Raw sequences for the cyanobacteria strain 193
R1=/scratch/project_2005590/COURSE_FILES/RAWDATA_ILLUMINA/Oscillatoriales-193_1.fastq.gz
R2=/scratch/project_2005590/COURSE_FILES/RAWDATA_ILLUMINA/Oscillatoriales-193_2.fastq.gz
```


You can check if your variable was set correctly by using:

```bash
echo $R1
echo $R2
```


### Running fastQC
Run `fastQC` to the files stored in the RAWDATA folder. What does the `-o` and `-t` flags refer to?

```bash
fastqc $R1 -o fastqc_raw/ -t 1

fastqc $R2 -o fastqc_raw/ -t 1
```

Copy the resulting HTML file to your local machine with `scp` from the command line (Mac/Linux) or *WinSCP* on Windows.  
Have a look at the QC report with your favourite browser.  

After inspecting the output, it should be clear that we need to do some trimming.  
__What kind of trimming do you think should be done?__

### Running Cutadapt

```bash
# To create a variable to your cyanobacterial strain:
strain=328
```

The adapter sequences that you want to trim are located after `-a` and `-A`.  
What is the difference with `-a` and `-A`?  
And what is specified with option `-p` or `-o`?
And how about `-m` and `-j`?  
You can find the answers from Cutadapt [manual](http://cutadapt.readthedocs.io).

Before running the script, we need to create the directory where the trimmed data will be written:

```bash
mkdir trimmed
```


```bash
cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -o trimmed/"$strain"_cut_1.fastq -p trimmed/"$strain"_cut_2.fastq $R1 $R2 --minimum-length 80 > cutadapt.log

```

### Running fastQC on the trimmed reads
You could now check the `cutadapt.log` and answer:

* How many read pairs we had originally?
* How many reads contained adapters?
* How many read pairs were removed because they were too short?
* How many base calls were quality-trimmed?
* Overall, what is the percentage of base pairs that were kept?

Then make a new folder (`FASTQC`) for the QC files of the trimmed data and run fastQC and multiQC again as you did before trimming:

```bash
mkdir fastqc_out_trimmed
fastqc trimmed/*.fastq -o fastqc_out_trimmed/ -t 1
```

Copy the resulting HTML file to your local machine as earlier and look how well the trimming went.  
Did you find problems with the sequences? We can further proceed to quality control using Prinseq.

#### Running Prinseq

You could check the different parameters that can be used in prinseq:
http://prinseq.sourceforge.net/manual.html

```bash

module load prinseq

```

Run program using the previous trimmed reads:

```bash
prinseq-lite.pl \
-fastq trimmed/"$strain"_cut_1.fastq \
-fastq2 trimmed/"$strain"_cut_2.fastq \
-min_qual_mean 25 \
-trim_left 10 \
-trim_right 8 \
-trim_qual_right 36 \
-trim_qual_left 30 \
-min_len 80 \
-out_good trimmed/"$strain"_pseq -log prinseq.log
```

You can check the `prinseq.log` and run again FastQC on the Prinseq trimmed sequences and copy them to your computer. You can now compare the quality of these sequences with the raw and cutadapt trimmed sequences FastQC results. Did you find any difference?

```bash
cd trimmed/

fastqc "$strain"_pseq_*.fastq -o ../fastqc_out_trimmed/ -t 1

```

### Optional - To compare raw and trimmed sequences using multiqc


To combine all the reports .zip in a new `combined_fastqc` folder with multiQC:
```bash
mkdir combined_fastqc

cp fastqc_raw/*zip combined_fastqc/
cp fastqc_out_trimmed/*zip combined_fastqc/

```

MultiQC is not pre-installed to Puhti, so we have created a virtual environment that has it.

```bash
export PROJAPPL=/projappl/project_2005590
module purge
module load bioconda/3
source activate mbdp_genomics
cd combined_fastqc/
multiqc . --interactive
```

To leave the interactive node, type `exit`.  

You can copy the file `multiqc_report.html` to your computer and open it in a web browser. Can you see any difference among the raw and trimmed reads?


## Genome assembly with Spades
Now that you have good trimmed sequences, we can assemble the reads.
For assembling you will need more resources than the default.  
Allocate 8 cpus, 20000 Mb of memory (20G) and 4 hours.  
Remember also the accounting project, `project_2005590`.

```bash

# Remember to modify  this
sinteractive --account --time --mem --cores


# Deactivate the current virtual environment and reset the modules before loading Spades
source deactivate mbdp_genomics
module purge


# Activate program
module load gcc/9.1.0
module load spades/3.15.0


# Cyano strain and processed reads
strain=328
R1=trimmed/"$strain"_pseq_1.fastq
R2=trimmed/"$strain"_pseq_2.fastq

```

### Run Spades

We will use the trimmed Illumina and Nanopore sequences to assemble the cyanobacteria genomes. Check the commands used using `spades.py -h`

```bash
spades.py --nanopore nanopore.trimmed.fastq.gz -1 $R1 -2 $R2 -o spades_hybrid_out -t 8
```

If you have time, you can try different options for assembly. Read more from [here](https://cab.spbu.ru/files/release3.15.0/manual.html) and experiment.  
Remember to rename the output folder for your different experiments.

After you're done, remember to close the interactive connection and free the resources with `exit`.



## Eliminate contaminant contigs with Kaiju

Kaiju is no pre-installed to Puhti so we need to reset the modules and activate the virtual environment again.

```bash
export PROJAPPL=/projappl/project_2005590
module purge
module load bioconda/3
source activate mbdp_genomics
```

### Run kaiju batch script
To run Kaiju you can use the script in `/scratch/project_2005590/COURSE_FILES/run_kaiju.sh`. This script takes your assembly as input and will eliminate all sequences not classified as cyanobacteria, creating a new "clean" file. You could go through the script and look at https://docs.csc.fi/computing/running/creating-job-scripts-puhti/ and `kaiju -h` before run in your folder:


```bash
sbatch /scratch/project_2005590/COURSE_FILES/run_kaiju.sh -i spades_hybrid_out/scaffolds.fasta -o kaiju_out
```

You can check the status of your job with:  

```bash
squeue -l -u $USER
```

After the job has finished, you can see how much resources it actually used and how many billing units were consumed.

```bash
seff JOBID
```

**NOTE:** Change **JOBID** the the job id number you got when you submitted the script.


### Optional - Visualizing the taxonomic assignment of contigs with the Kaiju web server

You can check the taxonomic assignment of the assembled contigs with some visuals using the [Kaiju web server](https://kaiju.binf.ku.dk/server).

This web tool only accepts compressed FASTA files (or FASTQ), so we need to compress our assembled genome file using `gzip`. The parameters on the web server can be left as is.

```bash
gzip -c your_assembly.fasta > your_assembly.fasta.gz
```


## Assembly QC

After the assembly has finished we will use Quality Assessment Tool for Genome Assemblies, [Quast](http://quast.sourceforge.net/) for (comparing and) evaluating our assemblies. Quast can be found from Puhti, but since there might be some incompability issues with Python2 and Python3, we will use a Singularity container that has Quast installed.  
More about Singularity: [More general introduction](https://sylabs.io/guides/3.5/user-guide/introduction.html) and [a bit more CSC specific](https://docs.csc.fi/computing/containers/run-existing/).

```
singularity exec --bind $PWD:$PWD /projappl/project_2005590/containers/quast_5.0.2.sif \
                                                  quast.py -o quast_out kaiju_out_filtered.fasta -t 4
```

Now you can move the file `quast_out/report.html` to your computer and look for the quality control files in the web browser of your preference.


## Calculate the genome coverage

To calculate the genome coverage, all the reads used for the assembly must be mapped to the final genome. For that, we can use three programs: Bowtie2 to map the reads; Samtools to sort and make an index of the mapped reads; and bedtools to make the calculation.

The entire workflow can take a long time, so the bedtools output with the sequencing depth for each base in the genome is available for each cyanobacterial strain in `/scratch/project_2005590/COURSE_FILES/RESULTS/coverage_[strain_number]/CoverageTotal.bedgraph` (don't forget to put the strain number).

You can check the commands used in the workflow to generate this file in the script in: `/scratch/project_2005590/COURSE_FILES/SCRIPTS/genome_coverage_workflow.sh`.


You can visualize the contents of the file `CoverageTotal.bedgraph` using the command `head` to show the first few lines.

The first 10 lines of `CoverageTotal.bedgraph` for the strain 328 as an example:
```bash
NODE_2_length_2022818_cov_20.647969   0   1   19
NODE_2_length_2022818_cov_20.647969   1   2   31
NODE_2_length_2022818_cov_20.647969   2   3   53
NODE_2_length_2022818_cov_20.647969   3   4   57
NODE_2_length_2022818_cov_20.647969   4   6   60
NODE_2_length_2022818_cov_20.647969   6   8   62
NODE_2_length_2022818_cov_20.647969   8   9   69
NODE_2_length_2022818_cov_20.647969   9   11  73
NODE_2_length_2022818_cov_20.647969   11  12  74
NODE_2_length_2022818_cov_20.647969   12  13  76
```
The first column refers to the contig name, the second and third column refers to the base position, and the fourth column is the sequencing depth of the base.

In order to calculate the final genome coverage we must calculate the average sequencing depth of all the bases in the genome. We can achieve this by using `awk`.

This command will sum all the numbers in the fourth column of the file `CoverageTotal.bedgraph` and divide by the total numbers of lines (number of bases), giving the average number. The final result will be printed in your screen (or you can save the result in a file using `> coverage.txt` in the end of the command).

```bash
cat CoverageTotal.bedgraph | awk '{total+=$4} END {print total/NR}'
```

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

## Genome annotation with Prokka

Now we can annotate our genome assembly using [Prokka](https://github.com/tseemann/prokka)

```bash
module purge
module load biokit
module load bioperl

prokka --cpus 8 --outdir prokka_out --prefix your_strain_name path-to/your_assembly.fasta
```

Check the files inside the output folder. Can you find the genes involved in the synthesis of Geosmin in one or more of these files?

### Optional - Annotation and visualization of CRISPR-Cas and Phages

Some genome features are better annotated when considering the genome context of a region involving many genes, instead of looking at only one gene at the time. Two examples of this case are the CRISPR-Cas system and Phages.

The CRISPR-Cas can be annotated using [CRISPRone](https://omics.informatics.indiana.edu/CRISPRone/denovo.php) and Phages can be annotated using [PHASTER](https://phaster.ca/).

Can you find any differences in the annotation of some specific genes when comparing the results of these tools with the Prokka annotation?


## Metaphlan

Next we will also analyze individual reads in addition to the assembly based approaches. Which files would you use for this? 

We will use a tool called [Metaphlan4](https://github.com/biobakery/biobakery/wiki/metaphlan4) to analyze these reads. 

# run gtdbtk
export GTDBTK_DATA_PATH=/scratch/project_2005590/databases/GTDB/release202/

singularity exec --bind $GTDBTK_DATA_PATH:$GTDBTK_DATA_PATH,$PWD:$PWD,$TMPDIR:/tmp  \
                    /projappl/project_2005590/containers/gtdbtk_1.7.0.sif \
                    gtdbtk classify_wf -x fasta --genome_dir PATH/TO/GENOME/FOLDER \
                    --out_dir OUTPUT/FOLDER --cpus 4 --tmpdir gtdb_test
```

## Pangenomics with Anvi'o


We will run the whole anvi'o part interactively. So again, allocate a computing node with enough memory (>40G) and 8 threads.  

Make a new folder for pangenomics
```
mkdir pangenomics
cd pangenomics
```

Make sure you have at least a copy of each of your genomes in one folder and make a new environmental variable pointing there.
```
GENOME_DIR=ABSOLUTE/PATH/TO/GENOME/DIR
CONTAINERS=/projappl/project_2005590/containers/
```

Then we will polish the contig names for each of the genomes before doing anything with the genomes.  
And also copy few annotated reference genomes to be included in our pangenome.  
```
singularity exec --bind $PWD:$PWD,$GENOME_DIR:/genomes $CONTAINERS/anvio_7.sif \
          anvi-script-reformat-fasta --simplify-names -o Oscillatoriales_193.fasta \
          -r reformat_193_report.txt /genomes/GENOME1.fasta

singularity exec --bind $PWD:$PWD,$GENOME_DIR:/genomes $CONTAINERS/anvio_7.sif \
          anvi-script-reformat-fasta --simplify-names -o Oscillatoriales_327_2.fasta \
          -r reformat_327_2_report.txt /genomes/GENOME2.fasta

singularity exec --bind $PWD:$PWD,$GENOME_DIR:/genomes $CONTAINERS/anvio_7.sif \
          anvi-script-reformat-fasta --simplify-names -o Oscillatoriales_328.fasta \
          -r reformat_328_report.txt /genomes/GENOME3.fasta

# copy reference genomes to your own folder
mkdir reference_genomes
cp /scratch/project_2005590/COURSE_FILES/closest_oscillatoriales_genomes/*.gbf ./reference_genomes
```

Although you already annotated your genomes, we'll do it once more, because we changed the names of the contigs.

```
module load biokit
for strain in $(ls *.fasta); do prokka --cpus 8 --outdir ./${strain%.fasta}_PROKKA --prefix ${strain%.fasta} $strain; done
```

Now we have both Genbank files from the reference genomes and also from our own genomes.  
Next things is to get the contigs, gene calls and annotations to separate files that anvi'o understands.

```
for genome in $(ls */*.gbf)
do
    singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif \
                                        anvi-script-process-genbank \
                                            -i $genome -O ${genome%.gbf} \
                                            --annotation-source prodigal \
                                            --annotation-version v2.6.3
done
```

The pangenomcis part will be done using a anvi'o workflow (read more from here: )  
And for that we need a file that specifies where the files from previous step are (called `fasta.txt`).  

```
echo -e "name\tpath\texternal_gene_calls\tgene_functional_annotation" > fasta.txt
for strain in $(ls */*-contigs.fa)
do
    strain_name=${strain#*/}
    echo -e ${strain_name%-contigs.fa}"\t"$strain"\t"${strain%-contigs.fa}"-external-gene-calls.txt\t"${strain%-contigs.fa}"-external-functions.txt"
done >> fasta.txt
```
In addition to the `fasta.txt` file we need also a configuration file.  
So make a file called `config.json` containing the following.  

```
{
    "workflow_name": "pangenomics",
    "config_version": "2",
    "max_threads": "8",
    "project_name": "Oscillatoriales_pangenome",
    "external_genomes": "external-genomes.txt",
    "fasta_txt": "fasta.txt",
    "anvi_gen_contigs_database": {
        "--project-name": "{group}",
        "--description": "",
        "--skip-gene-calling": "",
        "--ignore-internal-stop-codons": true,
        "--skip-mindful-splitting": "",
        "--contigs-fasta": "",
        "--split-length": "",
        "--kmer-size": "",
        "--skip-predict-frame": "",
        "--prodigal-translation-table": "",
        "threads": ""
    },
    "anvi_pan_genome": {
      "threads": "8"
    }
}
```
And then we're ready to run the whole pangenomics workflow.  

```
singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif anvi-run-workflow -w pangenomics -c config.json
```

When the workflow is ready, we can visualise the results interactively in anvi'o.  
For that we need to connect to Puhti bit differently.  
Log out from the current computing node, open a screen session and a new interactive connection to computing node.
Then take note of the computing node name, the login node number and your port number, you'll need them all in the next part.

### Tunneling the interactive interafce
Although you can install anvi'o on your own computer (and you're free to do so, but we won't have time to help in that), we will run anvi'o in Puhti and tunnel the interactive interface to your local computer.
To be able to to do this, everyone needs to use a different port for tunneling and your port number will be given on the course.

Connecting using a tunnel is a bit tricky and involves several steps, so pay special attention.
Detach from your screen and note on which login node you're on. Then re-attach and note the ID of the computing node your logged in. Then you will also need to remember your port number.

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

```
cd 03_PAN
export ANVIOPORT=PORT
singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif \
                                    anvi-display-pan \
                                        -g Oscillatoriales_pangenome-GENOMES.db \
                                        -p Oscillatoriales_pangenome-PAN.db \
                                        --server-only -P $ANVIOPORT

singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif \
                                    anvi-get-sequences-for-gene-clusters \
                                        -p Oscillatoriales_pangenome-PAN.db \
                                        -g Oscillatoriales_pangenome-GENOMES.db \
                                        -C default -b SCG \
                                        --concatenate-gene-clusters \
                                        -o single-copy-core-genes.fa                               

singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif \
                                    anvi-gen-phylogenomic-tree \
                                        -f single-copy-core-genes.fa  \
                                        -o SCG.tre

# study the geosmin phylogeny
singularity exec --bind $PWD:$PWD $CONTAINERS/anvio_7.sif \
                                    anvi-get-sequences-for-gene-clusters \
                                    -p Oscillatoriales_pangenome-PAN.db \
                                    -g Oscillatoriales_pangenome-GENOMES.db \
                                    --report-DNA-sequences \
                                    -C default -b Geosmin -o geosmin.fasta

# This time you need  to align the sequences (use  mafft) and construct the tree from  the aligned sequence file (raxml)
# concatenate the header before alignment
sed -i 's/[|:]/_/g' geosmin.fasta

module load biokit
ginsi geosmin.fasta > geosmin_aln.fasta

module purge
module load raxml
raxmlHPC -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s geosmin_aln.fasta -n geosmin
```

## Detection of secondary metabolites biosynthesis gene clusters

Biosynthetic genes putatively involved in the synthesis of secondary metabolites can identified using `antiSMASH`

Got to `https://antismash.secondarymetabolites.org/#!/start`. You can load the assembled genome you obtained and turn on all the extra features.

When the analysis is ready, you may be able to answer the following questions:

1. How many secondary metabolites biosynthetic gene clusters (BGC) were detected?
2. Which different types of BGC were detected and what is the difference among these types?
3. Do you think your strain produce all these metabolites? Why?


## Comparison of secondary metabolites biosynthesis gene clusters

Biosynthetic genes clusters can be compared using `BiG-SCAPE` (Biosynthetic Gene Similarity Clustering and Prospecting Engine) and `CORASON` (CORe Analysis of Syntenic Orthologs to prioritize Natural Product-Biosynthetic Gene Cluster)  https://bigscape-corason.secondarymetabolites.org/tutorial/index.html


You can run BiG-SCAPE/CORASON in your folder, using the results obtained from antiSMASH. The .gbk files from the three studied strains and selected reference strains from NCBI are available here: `/scratch/project_2005590/COURSE_FILES/BIGSCAPE/combinedGBK`

You need to load bigscape using a Conda environment:

```bash
export PROJAPPL=/projappl/project_2005590
module load bioconda/3
source activate bigscape
```


```bash
sinteractive -A project_2005590
```

Now you can run BiG-SCAPE/CORASON in your user folder:

```bash
mkdir bigscape
cd bigscape/

python /scratch/project_2005590/COURSE_FILES/BIGSCAPE/BiG-SCAPE/bigscape.py -c 4 -i /scratch/project_2005590/COURSE_FILES/BIGSCAPE/combinedGBK --pfam_dir /scratch/project_2005590/COURSE_FILES/BIGSCAPE/Pfam_database -o bigscape_auto --anchorfile /scratch/project_2005590/COURSE_FILES/BIGSCAPE/anchor_domains.txt --mode auto --hybrids-off --mibig --cutoffs 0.40 0.50 0.60
```

Take a look what it means each parameter used: https://git.wageningenur.nl/medema-group/BiG-SCAPE/-/wikis/home

Once the run is finished, you may transfer the folder to your own computer and observe the results.

Can you find the geosmin and/or 2-methylisoborneol biosynthetic gene clusters?




## Sandbox
Place to store some scratch code while testing.

### checkM
checkM should work from singularity container. Need to pull the right container (tag: 1.1.3--py_0) to course folder and test it once again
```
# needs computing node, otherwise runs out of memory (40G)
singularity exec --bind checkM_test/:/checkM_test /projappl/project_2005590/containers/checkM_1.1.3.sif \
              checkm lineage_wf -x fasta /checkM_test /checkM_test -t 4 --tmpdir /checkM_test
```
### GTDB-tk
Download database before running, needs >200G
```
# download gtdb database
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar -xzf gtdbtk_data.tar.gz

# run gtdbtk
export GTDBTK_DATA_PATH=/scratch/project_2005590 /databases/GTDB/release202/
singularity exec --bind $GTDBTK_DATA_PATH:$GTDBTK_DATA_PATH,$PWD:$PWD  /projappl/project_2005590/containers/gtdbtk_1.7.0.sif \
              gtdbtk classify_wf -x fasta --genome_dir checkM_test/ --out_dir gtdb_test --cpus 4  --tmpdir gtdb_test
```

### basecalling
```

~/projappl/ont-guppy/bin/guppy_basecaller -i fast5_pass/ -s BASECALLED/ -c ~/projappl/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg --device auto --min_qscore 10
cat BASECALLED/pass/*.fastq |gzip > BASECALLED/strain_328_nanopore.fastq.gz
```
