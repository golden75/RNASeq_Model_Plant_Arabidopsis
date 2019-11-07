# RNA-Seq for Model Plant (Arabidopsis thaliana)

This repository is a usable, publicly available tutorial for analyzing differential expression data and creating topological gene networks. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bio Informatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), and [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.   

  


Contents  
1.  [Introduction and programs]()   
2.  [Accessing the data using sra-toolkit]()   
3.  [Quality control using sickle]()   
4.  [Aligning reads to a genome using hisat2]()   
5.  [Transcript assembly and quantification with StringTie]()   
6.  [Differential expression analysis using ballgown]()   
7.  [Topological networking using cytoscape]()   
8.  [Conclusion]()   
  
  
## 1. Introduction and programs   

In this tutorial, we will be analyzing thale cress (Arabidopsis thaliana) RNA-Seq data from various parts of the plant (roots, stems). Perhaps one of the most common organisms for genetic study, the aggregrate wealth of genetic  Information of the thale cress makes it ideal for new-comers to learn. Organisms such as this we call "model organisms". You may think of model organisms as a subset of living things which, under the normal conventions of analysis, behave nicely. The data we will be analyzing comes from an experiment in which various cellular RNA was collected from the roots and shoots of a single thale cress. The RNA profiles are archived in the SRA, and meta Information on each may be viewed through the SRA ID: [SRR3498212](https://www.ncbi.nlm.nih.gov/sra?term=SRX1756762), [SRR3498213](https://www.ncbi.nlm.nih.gov/sra/?term=SRR3498213), [SRR3498215](https://www.ncbi.nlm.nih.gov/sra?term=SRX1756765), [SRR3498216](https://www.ncbi.nlm.nih.gov/sra?term=SRX1756766).    


The Single Read Archive, or SRA, is a publicly available database containing read sequences from a variety of experiments. Scientists who would like their read sequences present on the SRA submit a report containing the read sequences, experimental details, and any other accessory meta-data.

Our data, SRR3498212, SRR3498213, SRR3498215, SRR3498216 come from root 1, root 2, shoot 1, and shoot 2 of a single thale cress, respectively. Our objective in this analysis is to determine which genes are expressed in all samples, quantify the expression of each common gene in each sample, identify genes which are lowly expressed in roots 1 and 2 but highly expressed in shoots 1 and 2, or vice versa, quantify the relative expression of such genes, and lastly to create a visual topological network of genes with similar expression profiles.

You may connect to Xanadu via SSH, which will place you in your home directory   

```bash
cd /home/CAM/$USER 
```   

Your home directory contains adequate storage and will not pollute the capacities of other users on the cluster. To check the storage limits of your some directory please refer the [Xanadu](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/xanadu/) guide.   

The workflow may be cloned into the appropriate directory using the terminal command:  
```bash
git clone <git-hub-link>  
```  

Which will produce a copy of the working directory in your work space:  
```
RNASeq_Model_Plant_Arabidopsis
├── 01_raw_reads
└── tmp
```   

All of the completed scripts for this tutorial are available for you to submit in each directory. However, before submitting, you may want to edit the scripts to include your email!  

The tutorial will be using SLURM schedular to submit jobs to Xanadu cluster. In each script we will be using it will contain a header section which will allocate the resources for the SLURM schedular. The header section will contain:   

```bash
#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err  
```   

Before begining you need to create a directory called `tmp` in your home directory, which can be done using the following command:  
```bash
mkdir $HOME/tmp/
```  

As a precautionary measure, always include your temporary directory in the environment, using `export` command in each of your scripts. 
```bash
export TMPDIR=$HOME/tmp
```  

While not all programs require a temporary directory to work, it takes far less time including ours in the environment than it is waiting for an error! 

Before beginning, we need to understand a few aspects of the Xanadu server. When first logging into Xanadu from your local terminal, you will be connected to the submit node. The submit node is the interface with which users on Xanadu may `submit` their processes to the desired compute nodes, which will run the process. Never, under any circumstance, run processes directly in the submit node. Your process will be killed and all of your work lost! This tutorial will not teach you shell script configuration to submit your tasks on Xanadu. Therefore, before moving on, read and master the topics covered in the [Xanadu Tutorial](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/xanadu/).   

Now that we have covered the introduction and objective of our analysis, we may begin!   
  
  
## 2. Accessing the data using sra-toolkit   

We know that the SRA contain the read sequences and accessory meta Information from experiments. Rather than downloading experimental data through a browser, we may use the [sratoolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/), "fastq-dump" function to directly dump raw read data into the current terminal directory. Let's have a look at this function (it is expected that you have read the Xanadu tutorial, and are familiar with loading modules):   

We will be working in the **01_raw_reads/** folder.  

``` 
module load sratoolkit/2.8.2 
fastq-dump SRR3498212
fastq-dump SRR3498213
fastq-dump SRR3498215
fastq-dump SRR3498216
```   

For our needs, we will simply be using the accession numbers to dump our experimental data into our directory. We know our accession numbers, so let's write a shell script to retrieve our raw reads. There are a variety of text editors available on Xanadu. My preferred text editor is "nano". Therefore, we will be using nano or vi to write our shell script.   

The full slurm script we have prepared is called, [data_dump.sh](/01_raw_reads/data_dump.sh). Now that we have our script saved, we submit it to the compute nodes with the following command:
```bash
sbatch data_dump.sh
```   

Now we wait until we receive an email that our process has finished. Once the download is done, you should be having the following fastq files in the directory as follows:  
```
01_raw_reads/
├── SRR3498212.fastq
├── SRR3498213.fastq
├── SRR3498215.fastq
└── SRR3498216.fastq
```  


Let's take a look at one of our files using `head` command:  
```bash
bash-4.2$ head SRR3498212.fastq

@SRR3498212.1 SN638:767:HC555BCXX:1:1108:2396:1996 length=50
NTCAATCGGTCAGAGCACCGCCCTGTCAAGGCGGAAGCAGATCGGAAGAG
+SRR3498212.1 SN638:767:HC555BCXX:1:1108:2396:1996 length=50
#&#60;DDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SRR3498212.2 SN638:767:HC555BCXX:1:1108:2934:1998 length=50
NCAGTTTTGACCAGATAGGTCTCGCTAAGAAGATTGAGAAGATCGGAAGA
+SRR3498212.2 SN638:767:HC555BCXX:1:1108:2934:1998 length=50
#&#60;&#60;&#60;BDFFEEHIIIHIHEHEHIIEHHHH?&#60;GCCCHCCGICHH1&#60;GHHHIC
@SRR3498212.3 SN638:767:HC555BCXX:1:1108:3860:2000 length=50
NTCGCTTCGTAAGCGAAAGGCCGCGAGTTCGAAGATCGGAAGAGCACACG
```   

We see that for our first three runs we have  Information about the sampled read including its length followed by the nucleotide read and then a "+" sign. The "+" sign marks the beginning of the corresponding scores for each nucleotide read for the nucleotide sequence preceding the "+" sign.  


## Quality control using sickle  

Sickle performs quality control on illumina paired-end and single-end short read data using a sliding window. As the window slides along the fastq file, the average score of all the reads contained in the window is calculated. Should the average window score fall beneath a set threshold, [sickle](https://github.com/najoshi/sickle/blob/master/README.md) determines the reads responsible and removes them from the run. After visiting the SRA pages for our data, we see that our data are single end reads. Let's find out what sickle can do with these:  

To check the `sickle` command and options you can load the module and try `sickle` in a terminal window:
```bash
module load sickle/1.33

sickle
```   

Which will give you the *sickle* command options:
```bash
Usage: sickle <command> [options]

Command
pe      paired-end sequence trimming
se      single-end sequence trimming

--help, display this help and exit
--version, output version  Information and exit
```  

Since we are having single-end reads, and to see the single-end options, type `sickle se` in a terminal window, and it will show the following options that the program have.  
```bash
Usage sickle se [options] -f <fastq sequence file> -t <quality type> -o <trimmed fastq file>

Options
-f, --fastq-file, Input fastq file (required)
-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)
-o, --output-file, Output trimmed fastq file (required)
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.
-x, --no-fiveprime, Don't do five prime trimming.
-n, --trunc-n, Truncate sequences at position of first N.
-g, --gzip-output, Output gzipped files.
--quiet, Don't print out any trimming  Information
--help, display this help and exit
--version, output version  Information and exit
```   

The quality may be any score from 0 to 40. The default of 20 is much too low for a robust analysis. We want to select only reads with a quality of 35 or better. Additionally, the desired length of each read is 50bp. Again, we see that a default of 20 is much too low for analysis confidence. We want to select only reads whose lengths exceed 45bp. Lastly, we must know the scoring type. While the quality type is not listed on the SRA pages, most SRA reads use the "sanger" quality type. Unless explicitly stated, try running sickle using the sanger qualities. If an error is returned, try illumina. If another error is returned, lastly try solexa.  

Lets put all of this together for our sickle script using our downloaded fastq files:

```bash
module load sickle/1.33

sickle se -f ../01_raw_reads/SRR3498212.fastq -t sanger -o trimmed_SRR3498212.fastq -q 35 -l 45
sickle se -f ../01_raw_reads/SRR3498213.fastq -t sanger -o trimmed_SRR3498213.fastq -q 35 -l 45
sickle se -f ../01_raw_reads/SRR3498215.fastq -t sanger -o trimmed_SRR3498215.fastq -q 35 -l 45
sickle se -f ../01_raw_reads/SRR3498216.fastq -t sanger -o trimmed_SRR3498216.fastq -q 35 -l 45
```

The full slurm script is called [sickle.sh](/02_trimmed_reads/sickle.sh), which can be found in the **02_trimmed_reads** folder. Once the files have been trimmed you can find the following file structure inside the directory.  
```
02_trimmed_reads/
├── trimmed_SRR3498212.fastq
├── trimmed_SRR3498213.fastq
├── trimmed_SRR3498215.fastq
└── trimmed_SRR3498216.fastq
```   

It is helpful to see how the quality of the data has changed after using sickle. To do this, we will be using the commandline versions of [fastqc](https://www.bio Informatics.babraham.ac.uk/projects/fastqc/INSTALL.txt) and [MultiQC](http://multiqc. Info/docs/), These two programs simply create reports of the average quality of our trimmed reads, with some graphs. There is no way to view a --help menu for these programs in the command-line. However, their use is quite simple, we simply run "fastqc <trimmed_fastq>" or "multiqc -f -n trimmed trimmed". Do not worry too much about the options for MultiQC! Let's write our script: 
 
