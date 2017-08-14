---
title: "RNA-Seq workflow"
author: "Mary Piper, Meeta Mistry, Bob Freeman"
date: "Wednesday, March 8, 2017"
---

Approximate time: 90 minutes

## Learning Objectives:

* Describe and implement the RNA-Seq workflow to align reads to the reference genome 
* Describe tools and methods within the RNAseq work-flow
* Assessing input and output filetypes

## Running a Workflow

### Setting up

To get started with this lesson, we will start an interactive session and ask for 6 cores, by adding `-n 6` to the bsub command:

```bash
$ bsub -Is -n 6 -q interactive bash	
```

Change directories into the `unix_workshop` directory and copy the `reference_data` folder into your project directory:

```bash
$ cd ~/unix_workshop/rnaseq
```

You should have a directory tree setup similar to that shown below. it is best practice to have all files you intend on using for your workflow present within the same directory. In our case, we have our original FASTQ files and post-trimming data generated in the previous section. We also have all reference data files that will be used in downstream analyses.

```
rnaseq/
	├── raw_data/
	├── reference_data/
	   └── chr1.fa
	   └── chr1-hg19_genes.gtf
	├── meta/
	├── results/
	├── scripts/
	└── logs/
```

We previously described a general overview of the steps involved in RNA-Seq analysis, and in this session we will take our clean reads and align them to the reference genome.

![workflow](../img/workflow_alignment.png)


So let's get started by loading up some of the modules for tools we need for this section: 

```bash
 $ module load seq/STAR/2.5.2b seq/samtools/1.3 seq/htseq/0.6.1p1
```

Create an output directory for our alignment files:

```bash
$ mkdir results/STAR
```

In the script, we will eventually loop over all of our files and have the cluster work on each one in serial, then in parallel. For now, we're going to work on just one to set up our workflow.  To start we will use the trimmed first replicate in the Mov10 overexpression group, `Mov10_oe_1.qualtrim25.minlen35.fq` 


**NOTE: if you did not follow the last section, please execute the following command:** (this will copy over the required files into your home directory.)

```bash
$ cp -r /groups/hbctraining/unix_workshop_other/trimmed_fastq data/

```

## Read Alignment
The alignment process consists of choosing an appropriate reference genome to map our reads against, and performing the read alignment using one of several splice-aware alignment tools such as [STAR](https://github.com/alexdobin/STAR) or [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) (HISAT2 is a successor to both HISAT and TopHat2). The choice of aligner is a personal preference and also dependent on the computational resources that are available to you.
 
For this workshop we will be using STAR (Spliced Transcripts Alignment to a Reference), an aligner designed to specifically address many of the challenges of RNAseq read mapping. STAR is shown to have **high accuracy** and outperforms other aligners by more than a **factor of 50 in mapping speed (but also requires quite a bit of memory**). 

### STAR Alignment Strategy

STAR is shown to have **high accuracy** and outperforms other aligners by more than a **factor of 50 in mapping speed(but also requires quite a bit of memory)**. The algorithm achieves this highly efficient mapping by performing a two-step process:

1. Seed searching
2. Clustering, stitching, and scoring

#### Seed searching

For every read that STAR aligns, STAR will search for the longest sequence that exactly matches the reference genome:

![STAR_step1](../img/alignment_STAR_step1.png)
	
The different parts of the read that are mapped separately are called 'seeds'. So the first MMP that is mapped to the genome is called *seed1*.

STAR will then search again for only the unmapped portion of the read to find the next longest sequence that exactly matches the reference genome, which will be *seed2*. 

![STAR_step2](../img/alignment_STAR_step2.png)

This sequential searching of only the unmapped portions of reads underlies the efficiency of the STAR algorithm. STAR uses an uncompressed suffix array (SA) to efficiently search for the longest matching portions of the read, this allows for quick searching against even the largest reference genomes. Other slower aligners use algorithms that often search for the entire read sequence before splitting reads and performing iterative rounds of mapping. More details on the algorithm itself can be found in the [STAR publication](http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635). 

**If STAR does not find an exact matching sequence** for each part of the read due to mismatches or indels, the seed will be extended.

![STAR_step3](../img/alignment_STAR_step3.png)

**If extension does not give a good alignment**, then the poor quality or adapter sequence (or other contaminating sequence) will be soft clipped.

![STAR_step4](../img/alignment_STAR_step4.png)


#### Clustering, stitching, and scoring

The separate seeds are stitched together to create a complete read by first clustering the seeds together based on proximity to a set of seeds that have good alignment scores and are not multi-mapping.

Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.). 

![STAR_step5](../img/alignment_STAR_step5.png)

## Running STAR

Aligning reads using STAR is a two step process:   

1. Create a genome index 
2. Map reads to the genome

> A quick note on shared databases for human and other commonly used model organisms. The Orchestra cluster has a designated directory at `/groups/shared_databases/` in which there are files that can be accessed by any user. These files contain, but are not limited to, genome indices for various tools, reference sequences, tool specific data, and data from public databases, such as NCBI and PDB. So when using a tool and requires a reference of sorts, it is worth taking a quick look here because chances are it's already been taken care of for you. 

```bash
$ ls -l /groups/shared_databases/igenome/
```

### Creating a genome index

For this workshop we are using reads that originate from a small subsection of chromosome 1 (~300,000 reads) and so we are using only chr1 as the reference genome. Therefore, we cannot use any of the ready-made indices available in the `/groups/shared_databases/` folder.

For this workshop, we have already indexed the reference genome for you as this can take a while. We have provided the code below that you would use to index the genome for your future reference, but please **do not run the code below**. For indexing the reference genome, a reference genome (FASTA) is required and an annotation file (GTF or GFF3) is suggested a more accurate alignment of the reads. 

The basic options to **generate genome indices** using STAR are as follows:

* `--runThreadN`: number of threads
* `--runMode`: genomeGenerate mode
* `--genomeDir`: /path/to/store/genome_indices
* `--genomeFastaFiles`: /path/to/FASTA_file (reference genome)
* `--sjdbGTFfile`: /path/to/GTF_file (gene annotation)
* `--sjdbOverhang`: readlength -1

```bash
** DO NOT RUN**
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles chr1.fa --sjdbGTFfile chr1-hg19_genes.gtf --sjdbOverhang 99
```

The basic options for **mapping reads** to the genome using STAR are as follows:

* `--runThreadN`: number of threads
* `--readFilesIn`: /path/to/FASTQ_file
* `--genomeDir`: /path/to/genome_indices
* `--outFileNamePrefix`: prefix for all output files

We will also be using some advanced options:
* `--outSAMtype`: output filetype (SAM default)
* `--outSAMUnmapped`: what to do with unmapped reads
* `--outSAMattributes`: SAM attributes

Note that default filtering is applied in which the maximum number of multiple alignments allowed for a read is set to 10. If a read exceeds this number there is no alignment output. To change the default you can use `--outFilterMultimapNmax`, but for this lesson we will leave it as default. The advanced parameters that we are going to use are described below:

More details on STAR and its functionality can be found in the [user manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf), we encourage you to peruse through to get familiar with all available options.

Now let's put it all together! The full STAR alignment command is provided below.

> If you like you can copy-paste it directly into your terminal. Alternatively, you can manually enter the command, but it is advisable to first type out the full command in a text editor (i.e. [Sublime Text](http://www.sublimetext.com/) or [Notepad++](https://notepad-plus-plus.org/)) on your local machine and then copy paste into the terminal. This will make it easier to catch typos and make appropriate changes. 

```bash
$ STAR --runThreadN 6 \
--genomeDir /groups/hbctraining/unix_workshop_other/reference_STAR \
--readFilesIn data/trimmed_fastq/Mov10_oe_1.qualtrim25.minlen35.fq \
--outFileNamePrefix results/STAR/Mov10_oe_1_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI NM MD AS
```

#### Exercise
How many files do you see in your output directory? Using the `less` command take a look at `Mov10_oe_1_Log.final.out` and answer the following questions:  

1. How many reads are uniquely mapped?
2. How many reads map to more than 10 locations on the genome?
3. How many reads are unmapped due to read length?


### SAM/BAM
The output we requested from STAR is a BAM file, and by default returns a file in SAM format. **BAM is a binary version of the SAM file, also known as Sequence Alignment Map format.** The SAM file is a tab-delimited text file that contains information for each individual read and its alignment to the genome. The file begins with an optional header (which starts with '@'), followed by an alignment section in which each line corresponds to alignment information for a single read. **Each alignment line has 11 mandatory fields** for essential mapping information and a variable number of fields for aligner specific information.

These fields are described briefly below, but for more detailed information the paper by [Heng Li et al](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) is a good start.

![SAM](../img/sam_bam.png)

![SAM](../img/sam_bam3.png)

Let's take a quick look at our alignment. To do so we first convert our BAM file into SAM format using samtools and then pipe it to the `less` command. This allows us to look at the contents without having to write it to file (since we don't need a SAM file for downstream analyses).

```bash
$ samtools view -h results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam | less

```
 
Scroll through the SAM file and see how the fields correspond to what we expected.

### Assess the alignment (visualization)

Index the BAM file for visualization with IGV:

```bash
$ samtools index results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam
```

Use _**FileZilla**_ to copy the following files to your local machine:
 
`results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam`

`results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam.bai` 

> **NOTE: You can also transfer files to your laptop using the command line**
>
> Similar to the `cp` command, there is a command that allows you to securely copy files between computers. The command is called `scp` and allows files to be copied to, from, or between different hosts. It uses ssh for data transfer and provides the same authentication and same level of security as `ssh`. 
>
> First, identify the location of the _origin file_ you intend to copy, followed by the _destination_ of that file. Since the original file is located on Orchestra, this requires you to provide remote host and login information.

```bash
$ scp user_name@transfer.orchestra.med.harvard.edu:/home/user_name/unix_workshop/rnaseq_project/results/Mov10_oe_1_Aligned.sortedByCoord.out.bam* /path/to/directory_on_laptop
```

**Visualize**

* Start [IGV](https://www.broadinstitute.org/software/igv/download) _You should have this previously installed on your laptop_
* Load the Human genome (hg19) into IGV using the dropdown menu at the top left of your screen. _Note: there is also an option to "Load Genomes from File..." under the "Genomes" pull-down menu - this is useful when working with non-model organisms_
* Load the .bam file using the **"Load from File..."** option under the **"File"** pull-down menu. *IGV requires the .bai file to be in the same location as the .bam file that is loaded into IGV, but there is no direct use for that file.*

![IGV screenshot](../img/igv_screenshot.png)

#### Exercise
Now that we have done this for one sample, let's try using the same commands to perform alignment on one of the control samples. Using `Irrel_kd_1_qualtrim25.minlen35.fq` walk through the alignment commands above. Copy over the resulting BAM and index file to your laptop and upload into IGV for visualization. 

1. How does the MOV10 gene look in the control sample in comparison to the overexpression sample?
2. Take a look at a few other genes by typing into the search bar. For example, PPM1J and PTPN22. How do these genes compare? 

### Counting reads
Once we have our reads aligned to the genome, the next step is to count how many reads have been mapped to each gene. Counting is done with a tool called [`htseq-count`](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html). The input files required for counting include the BAM file and an associated gene annotation file in GTF format. `htseq-count` works by **taking the alignment coordinates for each read and cross-referencing that to the coordinates for features described in the GTF**. Most commonly a feature is considered to be a gene, which is the union of all exons (which is a feature type) that map to that gene. There is no minimum overlap to determine whether or not a read is counted for a particular gene, rather it is the mode that the user chooses. 

There are three modes available and are listed below in order of stringency, with most conservative at the top:

1. intersection-strict
2. union
3. intersection non-empty


We will be using the 'union' mode as it is default and most commonly used. To find out more on the different modes and how they affect your output, take a look at the [manual](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)

Let's start by creating a directory for the output:

```bash
$ mkdir results/counts
```

In it's most basic form the `htseq` command requires only the BAM file and the GTF file. We will add in a few additional parameters including `--format` to indicate BAM file, and `--stranded reverse` to specify that we have a stranded library created via the dUTP method. By default htseq-count will _ignore any reads that map to multiple locations_ on the genome. This results in undercounting but also helps reduce false positives. While multi-mappers are a feature that cannot be modified, there is a parameter that allows the user to filter reads by specifying a minimum alignment quality. 

You will notice at the end of the command we have added a redirection symbol. Since htseq-count outputs results to screen, we need to re-direct it to file.

```bash
$ htseq-count --stranded reverse --format bam results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam data/reference_data/chr1-hg19_genes.gtf  >  results/counts/Mov10_oe_1.counts
```

***
#### Exercise
Take a look at the end of the file using the `tail` command. You should see a summary of how the reads were classified. 

1. How many reads were assigned as no_feature? Why would they be classified this way?
2. How many reads were found to map to multiple locations?

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson were derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
