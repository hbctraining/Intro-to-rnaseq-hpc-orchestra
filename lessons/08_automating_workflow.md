---
title: "Automating an RNA-Seq workflow"
author: "Bob Freeman, Meeta Mistry, Radhika Khetani"
date: "Tuesday, August 22, 2017"
---

## Learning Objectives:

* Automate a workflow by grouping a series of sequential commands into a script
* Modify and submit the workflow script to the cluster

### From Sequence reads to Count matrix

Once you have optimized all the tools and parameters using a single sample, you can write a script to run the analysis on all the samples at the same time.

This will ensure that you run every sample with the exact same parameters, and will enable you to keep track of all the tools and their versions. In addition, the script is like a lab notebook, in the future you (or your colleagues) can go back and check the workflow for methods, which enables efficiency and reproducibility.

#### Automating this Workflow with a bash script

The easiest way for you to be able to repeat this process (from running STAR through to getting counts) is to capture the steps that
you've performed for `Mov10_oe_1` in a bash script. And you've already learned how to do this in previous
lessons. 

***

**Exercise**

* Using the text editor on your laptop, start creating a script with the commands you used to do the FastQC and alignment. Also include the command you might use for running `featureCounts` on only 1 sample at a time. 

> Note: You can use your command history to retrieve the commands for each step. Don't forget the "shebang line". Make sure you include the `mkdir` commands you used to create new directories.

* If you wanted to submit this script as a job on Orchestra, what would you have to modify? 
* What would have to modify to run this script on a different file?

***

#### Granting our Workflow More Flexibility

Hopefully this exercise gives you a glimpse into of how you might get a bash script started for a workflow you have just finished optimizing. We need to make several changes to this script to make it efficient, so we are going to start with a clean slate (a new file).

##### More variables

**The first major change is allowing flexibility in the input fastq file.** To be able to provide an input to the shell script, we need to use so-called **Positional Parameters**. 

The variables $1, $2, $3,...$9 and so on are positional parameters in the context of a shell script command. Let's take a closer look with an example.

E.g. We can refer to the components of the following command `sh run_rnaseq.sh input.fq input.gtf 12` as numbered variables inside the `run_rnaseq.sh` script, as follows:

`$0` = run_rnaseq.sh
`$1` = input.fq
`$2` = input.gtf
`$3` = 12

*There can be virtually unlimited numbers of inputs to a shell script, but it is wise to only have a few inputs to avoid errors and confusion when running a script that used positional parameters.*

> [This is an example of a simple script that used the concept of positional parameters and the associated variables](http://steve-parker.org/sh/eg/var3.sh.txt). You should try this script out after the class to get a better handle on positional parameters for shell scripting.

Let's use this new concept we have learned in the brand new script we are writing. We want the first positional parameter to be the name of our fastq file, and we will talk a little later about how to make sure you specify this for the user. We could just use the `$1` through out the script when we want to refer to the fastq file, but this is not intuitive so we want to create a new variable called `fq` and copy the contents of `$1` into it.

```bash
fq=$1
```
> When we set up variables we do not use the `$` before it, but when we *use the variable*, we always have to have the `$` before it. E.g setting up the `fq` variable => `fq=$1`, but we will use it as `fastqc $fq`.

Next we'll initialize 2 more variables named "genome" and "gtf, these will contain the paths to where the reference files are stored. This makes it easier to modify the script for when you want to use a different genome, i.e. you'll just have to change the contents of these variable at the beginning of the script.
```
# directory with genome reference FASTA and index files + name of the gene annotation file
genome=/groups/hbctraining/unix_workshop_other/reference_STAR/
gtf=~/unix_workshop/rnaseq_project/data/reference_data/chr1-hg19_genes.gtf
```
Next, make sure you load all the modules for the script to run. This is important so your script can run independent of any "prep" steps that need to be run beforehand:
    
```
# set up the software environment...
module load seq/fastqc/0.11.3
module load seq/STAR/2.5.3a
module load seq/samtools/1.3
PATH=/opt/bcbio/centos/bin:$PATH 	# for using featureCounts if not already in $PATH
```

We'll keep the output directory creation, however, we will add the `-p` option this will make sure that `mkdir` will create the directory only if it does not exist, and it won't throw an error if it does exist.

```
# make all of the output directories
# The -p option means mkdir will create the whole path if it 
# does not exist and refrain from complaining if it does exist

mkdir -p ~/unix_workshop/rnaseq_project/results/fastqc/
mkdir -p ~/unix_workshop/rnaseq_project/results/STAR
mkdir -p ~/unix_workshop/rnaseq_project/results/counts
```

In the script, it is a good idea to use echo for debugging/reporting to the screen (you can also use `set -x`):
```
echo "Processing file $fq ..."
```
> `set -x` debugging tool will display the command being executed, before the results of the command. In case of an issue with the commands in the shell script, this type of debugging lets you quickly pinpoint the step that is throwing an error. Often, tools will display the error that caused the program to stop running, so keep this in mind for times when you are running into issues where this is not availble.
> The command to turn it off is `set +x`

We also need to extract the "base name" of the file.
```
# grab base of filename for future naming
base=$(basename $fq .subset.fq)
echo "Sample name is $base"
```

> #### Remember `basename`?
>
> 1. the `basename` command: this command takes a path or a name and trims away all the information before the last `\` and if you specify the string to clear away at the end, it will do that as well. In this case, if the variable `$fq` contains the path *"unix_workshop_other/trimmed_fastq/Mov10_oe_1.qualtrim25.minlen35.fq"*, `basename $fq .qualtrim25.minlen35.fq` will output "Mov10_oe_1".
> 2. to assign this value to the `base` variable, we place the `basename...` command in parentheses and put a `$` outside. This syntax is necessary for assigning the output of a command to a variable.

Since we've already created our output directories, we can now specify all of our
output files in their proper locations. We will assign various file names to
 variables both for convenience but also to make it easier to see what 
is going on in the command below.

```
# set up output filenames and locations
fastqc_out=~/unix_workshop/rnaseq_project/results/fastqc/
align_out=~/unix_workshop/rnaseq_project/results/STAR/${base}_
counts_input_bam=~/unix_workshop/rnaseq_project/results/STAR/${base}_Aligned.sortedByCoord.out.bam
counts=~/unix_workshop/rnaseq_project/results/counts/${base}_featurecounts.txt
```
Our variables are now staged. We now need to modify the series of commands starting with STAR through to counts (htseq-count)
to use these variables so that it will run the steps of the analytical workflow with more flexibility:

```
# Run FastQC and move output to the appropriate folder
fastqc $fq
mv *fastqc* $fastqc_out

# Run STAR
STAR --runThreadN 6 --genomeDir $genome --readFilesIn $fq --outFileNamePrefix $align_out --outFilterMultimapNmax 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI NM MD AS

# Create BAM index
samtools index $counts_input_bam

# Count mapped reads
featureCounts -T 6 -s 2 -a $gtf -o $counts $counts_input_bam
```

It is always nice to have comments at the top of a more complex script to make sure that when your future self, or a co-worker, uses it they know exactly how to run it and what the script will do. So for our script, we can have the following lines of comments right at the top after `#!/bin/bash/`:

```
# This script takes a fastq file of RNA-Seq data, runs FastQC and outputs a counts file for it.
# USAGE: sh rnaseq_analysis_on_allfiles.sh <name of fastq file>
```

To transfer the saved file to Orchestra, you can copy and paste the script as a new file using `nano`.

```
$ cd ~/unix_workshop/rnaseq_project/scripts/

$ nano rnaseq_analysis_on_input_file.sh
```

*Alternatively, you can save the script on your computer and transfer it to `~/unix_workshop/rnaseq_project/scripts/` using FileZilla.*

Once the script has been saved, make it executable before running it. This is good to do even if your script runs fine without it; it will help avoid any future problems, and will enable your future self to know that it's an executable shell script.

```
$ chmod u+rwx rnaseq_analysis_on_input_file.sh 

$ sh rnaseq_analysis_on_input_file.sh <name of fastq>
```

#### Running our script iteratively as a job submission to the LSF scheduler

The above script will run in an interactive session **one file at a time**. But the whole point of writing this script was to run it on all files at once. How do you think we can do this?

To run the above script "in serial" for all of the files on a worker node via the job scheduler, we can create a separate submission script that will need 2 components:

1. **LSF directives** at the **beginning** of the script. This is so that the scheduler knows what resources we need in order to run our job on the compute node(s).
2. a for loop that iterates through and runs the above script for all the fastq files.

Below is what this second script would look like \[DO NOT RUN THIS\]:

```
#!/bin/bash

#BSUB -q training		# Partition to submit to (comma separated)
#BSUB -n 6              	# Number of cores, since we are running the STAR and featureCounts commands with 6 threads
#BSUB -W 1:30           	# Runtime in D-HH:MM (or use minutes)
#BSUB -R "rusage[mem=4000]"    	# Memory in MB
#BSUB -J rnaseq_mov10          	# Job name
#BSUB -o %J.out       		# File to which standard output will be written
#BSUB -e %J.err       		# File to which standard error will be written

# this `for` loop, will take the fastq files as input and run the script for all of them one after the other. 
for fq in ~/unix_workshop/rnaseq_project/raw_data/*.fq
do
  echo "running analysis on $fq"
  rnaseq_analysis_on_input_file.sh $fq
done
```

**But we don't want to run our the analysis on our 6 samples one after the other!** We want to run them "in parallel" as 6 separate jobs. 

**Note:** If you did run the above script, or something like it, you might want to give it the `.lsf` extension so it is obvious that it is meant to submit jobs to the LSF scheduler. 

***
**Exercise**

How would you run the above script?

***

#### Parallelizing workflow for efficiency

To parallelize our analysis, we will still need to write a second script, but with a modified way to specify the LSF directives. We will still be using a `for` loop. Parallelization will save you a lot of time with real (large) datasets.

Let's make a new file called `rnaseq_analysis_on_allfiles-for_lsf.sh`. *Note that this is a regular shell script.*

```bash
$ nano rnaseq_analysis_on_allfiles_for-lsf.sh
```

This file will loop through the same files as in the previous script, but the command it submits will be the actual bsub command:

```bash
#! /bin/bash

for fq in ~/unix_workshop/raw_fastq/*.fq
do
  
  bsub -q training -n 6 -W 1:30 -R "rusage[mem=4000]" -J rnaseq_mov10 -o %J.out -e %J.err "sh rnaseq_analysis_on_input_file.sh $fq"
  
  sleep 1	# wait 1 second between each job submission
  
done
```

In the above `for` loop please note that after the bsub directives the `sh rnaseq_analysis_on_input_file.sh $fq` command is in quotes.

> **NOTE:** All job schedulers are similar, but not the same. Once you understand how one works, you can transition to another one without too much trouble. They all have their pros and cons that the system administrators for a given cluster or HPC environment have taken into consideration and picked one that fits the needs of their users best. 

What you should see on the output of your screen would be the jobIDs that are returned from the scheduler for each of the jobs that your script submitted.

You can use the `bjobs` command to check progress (though there is a lag of about 60 seconds between what is happening and what is reported).

Don't forget about the `bkill` command, should something go wrong and you need to cancel your jobs.

#### Generating a Count Matrix

The above scripts will both generate count files for individual samples. Hence, after the script has run, you will have to do some cleanup using the `paste` and `cut` commands to generate a full count matrix wherein the first column is gene names and the rest of the columns are gene counts in each sample (with a header/column name). 

<img src="../img/count_matrix.png" width=500>

Alternatively, you could remove featureCounts from the script, and instead run it on all the BAM files after the jobs finish generating them.

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
