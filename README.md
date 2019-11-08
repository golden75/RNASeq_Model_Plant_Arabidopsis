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


## 3. Quality control using sickle  

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

It is helpful to see how the quality of the data has changed after using sickle. To do this, we will be using the commandline versions of [fastqc](https://www.bioInformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt) and [MultiQC](http://multiqc.Info/docs/), These two programs simply create reports of the average quality of our trimmed reads, with some graphs. There is no way to view a --help menu for these programs in the command-line. However, their use is quite simple, we simply run `fastqc <trimmed_fastq>` or `multiqc -f -n trimmed trimmed`. Do not worry too much about the options for MultiQC! Let's write our script:  

```bash
module load fastqc/0.11.7
module load MultiQC/1.7

fastqc trimmed_SRR3498212.fastq
fastqc trimmed_SRR3498213.fastq
fastqc trimmed_SRR3498215.fastq
fastqc trimmed_SRR3498216.fastq

multiqc -f -n trimmed trimmed*
```   

The full slurm script is called [quality_check.sh](/02_trimmed_reads/quality_check.sh), which can be found in the **02_trimmed_reads/** folder. As now you know the script can be run using the `sbatch quality_check.sh` command as we previously discussed. 

fastqc will create the files "trimmed_file_fastqc.html". To have a look at one, we need to move all of our "trimmed_file_fastqc.html" files into a single directory, and then [secure copy](https://www.techrepublic.com/article/how-to-use-secure-copy-for-file-transfer/) that folder to your local directory. Then, we may open our files! If that seems like too much work for you, you may open the files directly through this github. Simply click on any "html" file and you may view it in your browser immediately. Because of this, the steps mentioned above will not be placed in this tutorial.

This script will also create a directory "trimmed_data", which will hold text files and log files.   

```bash
02_trimmed_reads/
├── quality_check_483431.err
├── quality_check_483431.out
├── quality_check.sh
├── trimmed_data
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc.log
│   └── multiqc_sources.txt
├── trimmed.html
├── trimmed_SRR3498212_fastqc.html
├── trimmed_SRR3498212_fastqc.zip
├── trimmed_SRR3498213_fastqc.html
├── trimmed_SRR3498213_fastqc.zip
├── trimmed_SRR3498215_fastqc.html
├── trimmed_SRR3498215_fastqc.zip
├── trimmed_SRR3498216_fastqc.html
└── trimmed_SRR3498216_fastqc.zip
```

Let's have a look at the file format from fastqc and multiqc. When loading the fastqc file, you will be greeted with this screen:  
![](/images/fastqc1.png)  

There are some basic statistics which are all pretty self-explanatory. Notice that none of our sequences fail the quality report! It would be concerning if we had even one because this report is from our trimmed sequence! The same thinking applies to our sequence length. Should the minimum of the sequence length be below 45, we would know that sickle had not run properly. Let's look at the next index in the file:  
![](/images/fastqc2.png)  

This screen is simply a [box-and-whiskers plot](https://en.wikipedia.org/wiki/Box_plot) of our quality scores per base pair. Note that there is a large variance and lower mean scores (but still about in our desired range) for base pairs 1-5. These are the primer sequences! I will leave it to you to ponder the behavior of this graph. If you're stumped, you may want to learn how [Illumina sequencing](https://www.illumina.com/techniques/sequencing.html) works.  

Our next index is the per sequence quality scores:  
![](/images/fastqc3.png)  


This index is simply the total number of base pairs (y-axis) which have a given quality score (x-axis). This plot is discontinuous and discrete, and should you calculate the [Riemann sum](https://en.wikipedia.org/wiki/Riemann_sum) the result is the total number of base pairs present across all reads.  

The last index at which we are going to look is the "Overrepresented Sequences" index: 
![](/images/fastqc4.png)  

This is simply a list of sequences which appear disproportionately in our reads file. The reads file actually includes the primer sequences for this exact reason. When fastqc calculates a sequence which appears many times beyond the expected distribution, it may check the primer sequences in the reads file to determine if the sequence is a primer. If the sequence is not a primer, the result will be returned as "No Hit". Sequences which are returned as "No Hit" are most likely highly expressed genes.

We see that our multiqc file has the same indices as our fastqc files, but is simply the mean of all the statistics across our fastqc files:
![](/images/multiqc.png)  

  
 
## 4. Aligning reads to a genome using hisat2  

[HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml) is a fast and sensitive aligner for mapping next generation sequencing reads against a reference genome. HISAT2 requires two arguments: the reads file being mapped and the indexed genome to which those reads are mapped. Typically, the hisat2-build command is used to make a HISAT index file for the genome. It will create a set of files with the suffix .ht2, these files together build the index. What is an index and why is it helpful? Genome indexing is the same as indexing a tome, like an encyclopedia. It is much easier to locate  Information in the vastness of an encyclopedia when you consult the index, which is ordered in an easily navigatable way with pointers to the location of the  Information you seek within the encylopedia. Genome indexing is thus the structuring of a genome such that it is ordered in an easily navigatable way with pointers to where we can find whichever gene is being aligned. Let's have a look at how the hisat2-build command works:  

```bash
hisat2-build

HISAT2 version 2.1.0 by Daehwan Kim ( Infphilo@gmail.com, http://www.ccb.jhu.edu/people/ Infphilo)
Usage: hisat2-build [options] <reference_in> <ht2_index_base>
        reference_in            comma-separated list of files with ref sequences
        hisat2_index_base       write ht2 data to files with this dir/basename

-p : number of processors been used
--dta: report alignments tailored for transcript assemblers
-x: path to index generated from previous step
-q: query input files in fastq format
-S: output SAM file
``` 

As you can see, we simply enter our reference genome files and the desired prefix for our .ht2 files. Now, fortunately for us, Xanadu has many indexed genomes which we may use. To see if there is a hisat2 **Arabidopsis thaliana** indexed genome we need to look at the [Xanadu databases](https://bio Informatics.uconn.edu/databases/) page. We see that our desired indexed genomme is in the location /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/. Now we are ready to align our reads using hisat2 (for hisat2, the script is going to be writtten first with an explanation of the options after).  

```bash
module load hisat2/2.1.0

hisat2 -p 16 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q ../02_trimmed_reads/trimmed_SRR3498212.fastq  -S athaliana_root_1.sam
hisat2 -p 6 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q ../02_trimmed_reads/trimmed_SRR3498213.fastq -S athaliana_root_2.sam
hisat2 -p 6 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q ../02_trimmed_reads/trimmed_SRR3498215.fastq -S athaliana_shoot_1.sam
hisat2 -p 6 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q ../02_trimmed_reads/trimmed_SRR3498216.fastq -S athaliana_shoot_2.sam
```   

The full slurm script is called [hisat2_run.sh](/03_align/hisat2_run.sh) can be found at **03_align/** directory.  

Once the mapping have been completed, the file structure is as follows:
```
03_align/
├── athaliana_root_1.sam
├── athaliana_root_2.sam
├── athaliana_shoot_1.sam
└── athaliana_shoot_2.sam
```

When HISAT2 completes its run, it will summarize each of it’s alignments, and it is written to the standard error file, which can be found in the same folder once the run is completed. 

```
34475799 reads; of these:
  34475799 (100.00%) were unpaired; of these:
    33017550 (95.77%) aligned 0 times
    1065637 (3.09%) aligned exactly 1 time
    392612 (1.14%) aligned >1 times
4.23% overall alignment rate
42033973 reads; of these:
  42033973 (100.00%) were unpaired; of these:
    40774230 (97.00%) aligned 0 times
    931377 (2.22%) aligned exactly 1 time
    328366 (0.78%) aligned >1 times
3.00% overall alignment rate
31671127 reads; of these:
  31671127 (100.00%) were unpaired; of these:
    31103167 (98.21%) aligned 0 times
    465131 (1.47%) aligned exactly 1 time
    102829 (0.32%) aligned >1 times
1.79% overall alignment rate
49890217 reads; of these:
  49890217 (100.00%) were unpaired; of these:
    48622480 (97.46%) aligned 0 times
    1029943 (2.06%) aligned exactly 1 time
    237794 (0.48%) aligned >1 times
2.54% overall alignment rate
```  

Let's have a look at a SAM file:
```bash
 head -n 20 rnaseq_athaliana_root_1.sam
```  

```
@HD     VN:1.0  SO:unsorted
@SQ     SN:Chr1 LN:30427671
@SQ     SN:Chr2 LN:19698289
@SQ     SN:Chr3 LN:23459830
@SQ     SN:Chr4 LN:18585056
@SQ     SN:Chr5 LN:26975502
@SQ     SN:ChrM LN:366924
@SQ     SN:ChrC LN:154478
@PG     ID:hisat2       PN:hisat2       VN:2.1.0        CL:"/isg/shared/apps/hisat2/2.1.0/hisat2-align-s --wrapper basic-0 -p 16 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q trimmed_SRR3498212.fastq -S rnaseq_athaliana_root_1.sam"
SRR3498212.6    4       *       0       0       *       *       0       0       TTTCCAAGCCCTTTCTAGTCTGCGCTTGAGTTTGATTGCAGAGATCGGAA      DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      YT:Z:UU
SRR3498212.1    4       *       0       0       *       *       0       0       CAATCGGTCAGAGCACCGCCCTGTCAAGGCGGAAGCAGATCGGAAGAG        DDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII        YT:Z:UU
SRR3498212.4    4       *       0       0       *       *       0       0       AAAGGGCGTGGGTTCAAATCCCACAGATGTCACCAGATCGGAAGAGC DDHIIIIIIIIIEHHHIHIIIIHIIIIIIIIIIIIIIIIIIIIIIHH YT:Z:UU
SRR3498212.8    4       *       0       0       *       *       0       0       TTAAGATTGCTGATTTTGGCCTGGCACGTGAGGTTAAGATCGGAAGAGCA      DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      YT:Z:UU
SRR3498212.19   4       *       0       0       *       *       0       0       TGGATGATGGAAAAACCAGCAAGCCCCTCTTCTTTCAAGATCGGAAGAGC      DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      YT:Z:UU
SRR3498212.23   4       *       0       0       *       *       0       0       TTTGCCTTCCAAGCAATAGACCCGGGTAGATCGGAAGAGCACACGTCTGA      DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      YT:Z:UU
SRR3498212.24   4       *       0       0       *       *       0       0       TGAAACTTCTTGGTTTTAAAGTGTGAATATAGCTGACAAAAGATTGGAAG      DDDDDIIIIIIIIIIIIIIIIIIIHIIIIIIIIIIIIIIIIIIIIIIIII      YT:Z:UU
SRR3498212.12   4       *       0       0       *       *       0       0       AAGGGTGTTCTCTGCTACGGACCTCCAGATCGGAAGAGCACACGTCTGAA      DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      YT:Z:UU
SRR3498212.27   4       *       0       0       *       *       0       0       ATTGTTCCGGGCTGCCCAGTCCAAGCTGAGAGTGAAGATCGGAAGAGCAC      DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      YT:Z:UU
SRR3498212.29   4       *       0       0       *       *       0       0       TATGTCTACGCTGGTTCAAATCCAGCTCGGCCCACCAAGATCGGAAGAGC      DDDDDIIIIIIIIIHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      YT:Z:UU
SRR3498212.18   4       *       0       0       *       *       0       0       CGTGGGTTCGACTCCCACTGTGGTCGCCAAGATCGGAAGAGCACACGTC       DDDCHCCHHHEIHIGIIIEGHHIIIIGHHHIIIIIIIIIIIIIIIIIII       YT:Z:UU
```   

All of the lines starting with an "@" symbol tell us something about the chromosomes or our input. For instance "@SQ SN:Chr1 LN:30427671" tells us that we have a sequence (@SQ) whose sequence name is Chr1 (SN:Chr1), lastly the sequence has a length of 30427671bp (LN:30427671). You may be wondering what the first line means. It is quite straightfoward! The first line is simply the header (@HD) stating that the file is unsorted (SO:unsorted). The second column in the first line is somewhat of a dummy variable, but stands for "version number". Lastly we have the "@PG" line, which, in order, keeps track of the software used to write the file (ID:hisat2), the program name used to align the reads (PN:hisat2), the version of the program used (VN:2.1.0), and lastly the user input which started the process (written in the form that the program reads, not in which we wrote it).

The alignment portion of the SAM file is much more straight-forward and may be understood by reading the SAM output formatting guide linked in the beginning of this tutorial.

Because of the density of the sam file, it is compressed to binary to create a more easily tractable file for manipulation by future programs. We convert the sam file to bam with the following command and sort it such that the alignments are listed in the order the genes appear in the genome. To do this we use the software [samtools](https://en.wikipedia.org/wiki/SAMtools):  

SAMtools options can be looked at as:  
```bash
module load samtools/1.9

Usage :   samtools <command> [options]

Commands
  -- Indexing
     dict           create a sequence dictionary file
     faidx          index/extract FASTA
     index          index alignment

  -- Editing
     calmd          recalculate MD/NM tags and '=' bases
     fixmate        fix mate  Information
     reheader       replace BAM header
     targetcut      cut fosmid regions (for fosmid pool only)
     addreplacerg   adds or replaces RG tags
     markdup        mark duplicates

  -- File operations
     collate        shuffle and group alignments by name
     cat            concatenate BAMs
     merge          merge sorted alignments
     mpileup        multi-way pileup
     sort           sort alignment file
     split          splits a file by read group
     quickcheck     quickly check if SAM/BAM/CRAM file appears intact
     fastq          converts a BAM to a FASTQ
     fasta          converts a BAM to a FASTA

  -- Statistics
     bedcov         read depth per BED region
     depth          compute the depth
     flagstat       simple stats
     idxstats       BAM index stats
     phase          phase heterozygotes
     stats          generate stats (former bamcheck)

  -- Viewing
     flags          explain BAM flags
     tview          text alignment viewer
     view           SAM<->BAM<->CRAM conversion
     depad          convert padded BAM to unpadded BAM
```  

We are truly only interested in sorting our SAM files.

```bash
samtools sort

Usage: samtools sort [options...] [in.bam]

Options:
  -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
  -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
  -n         Sort by read name
  -t TAG     Sort by value of TAG. Uses position as secondary index (or read name if -n is set)
  -o FILE    Write final output to FILE rather than standard output
  -T PREFIX  Write temporary files to PREFIX.nnnn.bam
      --input-fmt-option OPT[=VAL]
               Specify a single input file format option in the form
               of OPTION or OPTION=VALUE
  -O, --output-fmt FORMAT[,OPT[=VAL]]...
               Specify output format (SAM, BAM, CRAM)
      --output-fmt-option OPT[=VAL]
               Specify a single output file format option in the form
               of OPTION or OPTION=VALUE
      --reference FILE
               Reference sequence FASTA FILE [null]
  -@, --threads INT
               Number of additional threads to use [0]
```


The sort function converts SAM files to BAM automatically. Therefore, we can cut through most of these options and do a simple `samtools sort -o <output.bam> <inupt.sam>`. Let's write our script:  

```bash
module load samtools/1.9

samtools sort -@ 6 -o athaliana_root_1.bam athaliana_root_1.sam
samtools sort -@ 6 -o athaliana_root_2.bam athaliana_root_2.sam
samtools sort -@ 6 -o athaliana_shoot_1.bam athaliana_shoot_1.sam
samtools sort -@ 6 -o athaliana_shoot_2.bam athaliana_shoot_2.sam
```   
The full slurm script is called [samsort.sh](/03_align/samsort.sh) and its stored in **03_align/** directory. 

This will create sorted bam files and the file structure in the **03_align/** directory will look like:  
```
03_align/
├── athaliana_root_1.bam
├── athaliana_root_2.bam
├── athaliana_shoot_1.bam
└── athaliana_shoot_2.bam
```   


##  5. Transcript assembly and qunatification with StringTie   

This is the most intricate part of our analysis. Because of the many steps involved, this tutorial is going to walk through the steps before writing a final batch script to be submitted.

Analysis of RNA-seq experiments requires accurate reconstructions of all the [isoforms](https://en.wikipedia.org/wiki/Gene_isoform) expressed from each gene, as well as estimates of the relative abundance of those isoforms. Accurate quantification benefits from knowledge of precisely which reads originated from each isoform, which cannot be computed perfectly because reads are much shorter than transcripts. StringTie assembles transcripts from RNA-seq reads that have been aligned to the genome, first grouping the reads into distinct gene loci and then assembling each locus into as many isoforms as are needed to explain the data. To begin, it would be nice to have a source describing all genome features (simply a list of genes, exons, transcripts, etc., which states on which chromosome and over which base-pairs those genes are). The file format for which we're looking is GFF (which is exactly as described the sentence prior). We can download the GFF file for the thale cress with the following code, when executing the command to download, make sure you are using an interative session.

```bash
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
```  

Once downloaded lets have a look at the first few lines in the GFF file: 
```bash
head TAIR10_GFF3_genes.gff 

Chr1    TAIR10  chromosome      1       30427671        .       .       .       ID=Chr1;Name=Chr1
Chr1    TAIR10  gene    3631    5899    .       +       .       ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
Chr1    TAIR10  mRNA    3631    5899    .       +       .       ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1
Chr1    TAIR10  protein 3760    5630    .       +       .       ID=AT1G01010.1-Protein;Name=AT1G01010.1;Derives_from=AT1G01010.1
Chr1    TAIR10  exon    3631    3913    .       +       .       Parent=AT1G01010.1
Chr1    TAIR10  five_prime_UTR  3631    3759    .       +       .       Parent=AT1G01010.1
Chr1    TAIR10  CDS     3760    3913    .       +       0       Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1    TAIR10  exon    3996    4276    .       +       .       Parent=AT1G01010.1
Chr1    TAIR10  CDS     3996    4276    .       +       2       Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1    TAIR10  exon    4486    4605    .       +       .       Parent=AT1G01010.1
```  

The GFF file is quite self-explanatory. However, it'd be nice if could combine all of the pieces of  Information from the GFF into something better. For instance, if there are multiple overlapping, but distinct exons from a single gene, we could use that  Information to determine the isoforms of that gene. Then, we could make a file which gives each isoform its own track (there are other extrapolations to be made, but this is our most relevant example). Luckily for us, we can use the program "gffread" to transform our GFF file into the more useful form just stated, The output of [gffread --help](https://github.com/gpertea/gffread) is much too dense for us to go into here, but the necessary options will be explained. Do not run this code! We are compiling this code with various other chunks into one script, be patient!   


```bash
module load gffread/0.9.12
gffread TAIR10_GFF3_genes.gff -T -o athaliana_TAIR10_genes.gtf
```   

The option -T tells gffread to convert our input into the gtf format, and the option -o simply is how we call the output. The GTF format is simply the transcript assembly file, and is composed of exons and coding sequences. Let's have a look at the GTF file:  

```bash
head athaliana_TAIR10_genes.gtf
```   

```

Chr1    TAIR10  exon    3631    3913    .       +       .       transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1    TAIR10  exon    3996    4276    .       +       .       transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1    TAIR10  exon    4486    4605    .       +       .       transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1    TAIR10  exon    4706    5095    .       +       .       transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1    TAIR10  exon    5174    5326    .       +       .       transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1    TAIR10  exon    5439    5899    .       +       .       transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1    TAIR10  CDS     3760    3913    .       +       0       transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1    TAIR10  CDS     3996    4276    .       +       2       transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1    TAIR10  CDS     4486    4605    .       +       0       transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
Chr1    TAIR10  CDS     4706    5095    .       +       0       transcript_id "AT1G01010.1"; gene_id "AT1G01010"; gene_name "AT1G01010";
```  

To look at the last few lines 
```bash
tail athaliana_TAIR10_genes.gtf
```   

```
ChrM    TAIR10  exon    349830  351413  .       -       .       transcript_id "ATMG01360.1"; gene_id "ATMG01360"; gene_name "ATMG01360";
ChrM    TAIR10  CDS     349830  351413  .       -       0       transcript_id "ATMG01360.1"; gene_id "ATMG01360"; gene_name "ATMG01360";
ChrM    TAIR10  exon    360717  361052  .       -       .       transcript_id "ATMG01370.1"; gene_id "ATMG01370"; gene_name "ATMG01370";
ChrM    TAIR10  CDS     360717  361052  .       -       0       transcript_id "ATMG01370.1"; gene_id "ATMG01370"; gene_name "ATMG01370";
ChrM    TAIR10  exon    361062  361179  .       -       .       transcript_id "ATMG01380.1"; gene_id "ATMG01380"; gene_name "ATMG01380";
ChrM    TAIR10  exon    361350  363284  .       -       .       transcript_id "ATMG01390.1"; gene_id "ATMG01390"; gene_name "ATMG01390";
ChrM    TAIR10  exon    363725  364042  .       +       .       transcript_id "ATMG01400.1"; gene_id "ATMG01400"; gene_name "ATMG01400";
ChrM    TAIR10  CDS     363725  364042  .       +       0       transcript_id "ATMG01400.1"; gene_id "ATMG01400"; gene_name "ATMG01400";
ChrM    TAIR10  exon    366086  366700  .       -       .       transcript_id "ATMG01410.1"; gene_id "ATMG01410"; gene_name "ATMG01410";
ChrM    TAIR10  CDS     366086  366700  .       -       0       transcript_id "ATMG01410.1"; gene_id "ATMG01410"; gene_name "ATMG01410";
```   

We see that whereas in our GFF file we have various untranslated regions included, as well as annotations, the GTF format contains  Information only on various transcripts for each gene. The "transcript_id" denoter in the last column tells us the gene and its isoform, and everything else about the GTF file is quite apparent!

Just as was stated for our conversion from gff to gtf, it would be helpful for us to perform the same operation on our aligned reads. That is, if there are multiple, overlapping but distinct reads from a single gene, we could combine these reads into one transcript isoform. Because we have the gene isoforms in the gtf file, we can re-map each assembled transcript to a gene isoform and then count how many mappings there are per isoform. This, in effect, allows us to quantify the expression rates of each isoform. We will be using the program [StringTie](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) to assemble the transcripts for each sample. StringTie requires three input arguments: the BAM alignment file, the genomic GTF file, and the desired output GTF filename. Thus, our code will look like (do not run this!):

```bash
module load stringtie/2.0.3

stringtie -p 16 ../03_align/athaliana_root_1.bam -G athaliana_TAIR10_genes.gtf -o athaliana_root_1.gtf

stringtie -p 16 ../03_align/athaliana_root_2.bam -G athaliana_TAIR10_genes.gtf -o athaliana_root_2.gtf

stringtie -p 16 ../03_align/athaliana_shoot_1.bam -G athaliana_TAIR10_genes.gtf -o athaliana_shoot_1.gtf

stringtie -p 16 ../03_align/athaliana_shoot_2.bam -G athaliana_TAIR10_genes.gtf -o athaliana_shoot_2.gtf
```

Let's have a look at the stringtie output:

```bash
athaliana_root_1.gtf

# stringtie -p 16 rnaseq_athaliana_root_1.bam -G athaliana_TAIR10_genes.gtf -o rnaseq_athaliana_root_1.gtf
# StringTie version 1.3.3b
Chr1	StringTie	transcript	28500	28706	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; reference_id "AT1G01046.1"; ref_gene_id "AT1G01046"; ref_gene_name "AT1G01046"; cov "0.241546"; FPKM "3.727008"; TPM "0.747930";
Chr1	StringTie	exon	28500	28706	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; reference_id "AT1G01046.1"; ref_gene_id "AT1G01046"; ref_gene_name "AT1G01046"; cov "0.241546";
Chr1	StringTie	transcript	47494	48839	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; cov "3.928230"; FPKM "60.611832"; TPM "12.163484";
Chr1	StringTie	exon	47494	47982	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1"; cov "4.529652";
Chr1	StringTie	exon	48075	48839	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "2"; cov "3.543791";
Chr1	StringTie	transcript	50075	51199	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; reference_id "AT1G01100.2"; ref_gene_id "AT1G01100"; ref_gene_name "AT1G01100"; cov "8.437494"; FPKM "130.188904"; TPM "26.126097";
Chr1	StringTie	exon	50075	50337	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "1"; reference_id "AT1G01100.2"; ref_gene_id "AT1G01100"; ref_gene_name "AT1G01100"; cov "6.228601";
Chr1	StringTie	exon	50419	50631	1000	-	.	gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "2"; reference_id "AT1G01100.2"; ref_gene_id "AT1G01100"; ref_gene_name "AT1G01100"; cov "9.487524";
```   

While this is certainly confusing, we can still understand it. To start each row we have the chromosome for each sequence. The second column is the software used to assemble the transcript, in our case StringTie. The third column is the sequence type, either a transcript or an exon. The next two columns are the start and end bp of the feature (assuming the chromosome starts at bp 1), followed by another column which is the score. The column after the score is either "+" or "-" for forward strand and reverse strand, respectively. Our last two columns are the frame ("." means frame not determined) and the feature meta Information.

You may be wondering about the last three columns which are not standard in GTF. The column "cov" is simply the covariance of the gene across samples (if the gene is highly or lowly expressed in both samples the covariance will be high, if it is highly expressed in one sample and lowly expressed in another sample, the covariance will be low). The FPKM column is the fragment per kilobase million value. Simply put, this is the number of times the specific feature was counted divided by how many counts there were for all features combined, in millions. That number is then divided by the length of the feature in kilobases. The reason for this being that longer features should have higher counts. For instance, when we split our mRNA in sequences of 50 or less for reading, one 5000 bp transcript will appear as 100 counts, and one 500 bp transcript will appear as 10 counts. Now, let's divide each count by the transcript length in kilobases: we have 20 as the value for the 5000 bp sequence (100/(5000/1000)) and 20 as the value for the 500 bp sequence (10/(500/1000)! Now our quantification matches the actual expression profile -- that is, both features were transcribed the same amount.

The last column, TPM, is the transcripts per feature per  million transcripts counted across all features combined. As evident from the example above, without some sort of scaling factor, this value is highly misleading.

The assembled isoforms in each sample are most likely different from those of other samples. However, we may repeat the process of determining isoforms but this time using the gtf files for all four samples. That is, if there are isoforms across the sample gtfs which overlap but are distinct, we may merge those into a single isoform. To do this we will be using the --merge option of stringtie. The --merge option of stringtie requires three input arguments: the genomic GTF file, the desired output filename, and a plain-text file containing each file to be merged separated by a return. To begin, let's make our plain-text file:

```bash
mergelist.txt

athaliana_root_1.gtf
athaliana_root_2.gtf
athaliana_shoot_1.gtf
athaliana_shoot_2.gtf
```   

After saving our plain-text file, we can merge our samples using the following code (do not run this!):

```bash
module load stringtie/2.0.3 
stringtie --merge -p 16 -G -o merged stringtie_merged.gtf
```  

Then;
```bash
module load gffcompare/0.10.4
gffcompare -r athaliana_TAIR10_genes.gtf -o merge stringtie_merged.gtf
``` 

While the options are quite self-explanatory, one thing to note is that the "-o" option is simply the output *prefix*. After running, you should see the following files in your directory (but hopefully you listened and did not run it yet!):  

So the files which is produced:
```
merged.annotated.gtf
merged.loci
merged.stats
merged.stringtie_merged.gtf.refmap
merged.stringtie_merged.gtf.tmap
merged.tracking
```   

Although you have not run the code yet, let's have a look at some of the files we've generated (we will not be looking at the merged.annotated.gtf as we are already quite familiar with the gtf format!)

```bash
head merged.loci
```   

```
XLOC_000001     Chr1[+]3631-5899        AT1G01010|AT1G01010.1   AT1G01010.1
XLOC_000002     Chr1[+]23146-31227      AT1G01040|AT1G01040.1,AT1G01040|AT1G01040.2     AT1G01040.1,AT1G01040.2
XLOC_000003     Chr1[+]28500-28706      AT1G01046|AT1G01046.1   AT1G01046.1
XLOC_000004     Chr1[+]44677-44787      AT1G01073|AT1G01073.1   AT1G01073.1
XLOC_000005     Chr1[+]52239-54692      AT1G01110|AT1G01110.2,AT1G01110|AT1G01110.1     AT1G01110.2,AT1G01110.1
XLOC_000006     Chr1[+]56624-56740      AT1G01115|AT1G01115.1   AT1G01115.1
XLOC_000007     Chr1[+]72339-74096      AT1G01160|AT1G01160.1,AT1G01160|AT1G01160.2     AT1G01160.1,AT1G01160.2
XLOC_000008     Chr1[+]75583-76758      AT1G01180|AT1G01180.1   AT1G01180.1
XLOC_000009     Chr1[+]88898-89745      AT1G01210|AT1G01210.1   AT1G01210.1
XLOC_000010     Chr1[+]91376-95651      AT1G01220|AT1G01220.1   AT1G01220.1
```   

We see we have a condensed form of our various exons. The exon loci name is the first column, followed by the chromosome, strand, and bp location in the second column. The final columns are the gene ID, the transcript IDs to which the loci belong, and the isoforms to which the transcripts belong.  


**merged.stats**  

```
less merged.stats


# gffcompare v0.10.4 | Command line was:
# gffcompare -r athaliana_TAIR10_genes.gtf -G -o merged stringtie_merged.gtf
#

# Summary for dataset: stringtie_merged.gtf
#     Query mRNAs :   41854 in   33403 loci  (30272 multi-exon transcripts)
#            (6013 multi-transcript loci, ~1.3 transcripts per locus)
# Reference mRNAs :   41607 in   33350 loci  (30127 multi-exon)
# Super-loci w/ reference transcripts:    33184
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    99.9    |
        Exon level:   100.0     |    99.9    |
      Intron level:   100.0     |   100.0    |
Intron chain level:   100.0     |    99.5    |
  Transcript level:   100.0     |    99.4    |
       Locus level:   100.0     |    99.8    |

     Matching intron chains:   30127
       Matching transcripts:   41607
              Matching loci:   33350

          Missed exons:       0/169264  (  0.0%)
           Novel exons:      90/169578  (  0.1%)
        Missed introns:       0/127896  (  0.0%)
         Novel introns:      13/127951  (  0.0%)
           Missed loci:       0/33350   (  0.0%)
            Novel loci:      75/33403   (  0.2%)

 Total union super-loci across all input datasets: 33403
41854 out of 41854 consensus transcripts written in merged.annotated.gtf (0 discarded as redundant)
```   

The Information here is quite apparent.  

```bash
head merged.stringtie_merged.gtf.refmap
```   

```
ref_gene_id     ref_id  class_code      qry_id_list
AT1G01010       AT1G01010.1     =       AT1G01010|AT1G01010.1
AT1G01040       AT1G01040.1     =       AT1G01040|AT1G01040.1
AT1G01040       AT1G01040.2     =       AT1G01040|AT1G01040.2
AT1G01046       AT1G01046.1     =       AT1G01046|AT1G01046.1
AT1G01073       AT1G01073.1     =       AT1G01073|AT1G01073.1
AT1G01110       AT1G01110.2     =       AT1G01110|AT1G01110.2
AT1G01110       AT1G01110.1     =       AT1G01110|AT1G01110.1
AT1G01115       AT1G01115.1     =       AT1G01115|AT1G01115.1
AT1G01160       AT1G01160.1     =       AT1G01160|AT1G01160.1
```   

Here we have the gene IDs from the reference GFF, followed by the isoform IDs from the reference GTF, the "class_code" which simply tells you that the last column was matched fully ("=") or partially ("c"). Taking our first column, we see that all of isoform ATG01010.1 was matched to the gene ATG01010.

```bash
head merged.stringtie_merged.gtf.tmap
```

```
ref_gene_id     ref_id  class_code      qry_gene_id     qry_id  num_exons       FPKM    TPM             cov     len     major_iso_id    ref_match_len
AT1G01090       AT1G01090.1     =       AT1G01090       AT1G01090.1     3       0.000000        0.000000        0.000000        1627    AT1G01090.1     1627
AT1G01100       AT1G01100.2     =       AT1G01100       AT1G01100.2     4       0.000000        0.000000        0.000000        631     AT1G01100.2     631
AT1G01100       AT1G01100.1     =       AT1G01100       AT1G01100.1     4       0.000000        0.000000        0.000000        587     AT1G01100.2     587
AT1G01100       AT1G01100.4     =       AT1G01100       AT1G01100.4     4       0.000000        0.000000        0.000000        607     AT1G01100.2     607
AT1G01100       AT1G01100.3     =       AT1G01100       AT1G01100.3     5       0.000000        0.000000        0.000000        566     AT1G01100.2     566
AT1G01120       AT1G01120.1     =       AT1G01120       AT1G01120.1     1       0.000000        0.000000        0.000000        1899    AT1G01120.1     1899
AT1G01180       AT1G01180.1     =       AT1G01180       AT1G01180.1     1       0.000000        0.000000        0.000000        1176    AT1G01180.1     1176
AT1G01183       AT1G01183.1     =       AT1G01183       AT1G01183.1     1       0.000000        0.000000        0.000000        101     AT1G01183.1     101
AT1G01225       AT1G01225.1     =       AT1G01225       AT1G01225.1     2       0.000000        0.000000        0.000000        1025    AT1G01225.1     1025
```

All of the  Information in the .tmap file may be readily understood now, knowing that "len" and "ref_match_len" are the sequence lengths and reference lengths, respectively.  

Lastly,
```bash
head merged.tracking
```   

```
TCONS_00000001  XLOC_000001     AT1G01010|AT1G01010.1   =       q1:AT1G01010|AT1G01010.1|6|0.000000|0.000000|0.000000|1688
TCONS_00000002  XLOC_000002     AT1G01040|AT1G01040.1   =       q1:AT1G01040|AT1G01040.1|20|0.000000|0.000000|0.000000|6251
TCONS_00000003  XLOC_000002     AT1G01040|AT1G01040.2   =       q1:AT1G01040|AT1G01040.2|20|0.000000|0.000000|0.000000|5877
TCONS_00000004  XLOC_000003     AT1G01046|AT1G01046.1   =       q1:AT1G01046|AT1G01046.1|1|0.000000|0.000000|0.000000|207
TCONS_00000005  XLOC_000004     AT1G01073|AT1G01073.1   =       q1:AT1G01073|AT1G01073.1|1|0.000000|0.000000|0.000000|111
TCONS_00000006  XLOC_000005     AT1G01110|AT1G01110.2   =       q1:AT1G01110|AT1G01110.2|5|0.000000|0.000000|0.000000|1782
TCONS_00000007  XLOC_000005     AT1G01110|AT1G01110.1   =       q1:AT1G01110|AT1G01110.1|3|0.000000|0.000000|0.000000|1439
TCONS_00000008  XLOC_000006     AT1G01115|AT1G01115.1   =       q1:AT1G01115|AT1G01115.1|1|0.000000|0.000000|0.000000|117
TCONS_00000009  XLOC_000007     AT1G01160|AT1G01160.1   =       q1:AT1G01160|AT1G01160.1|5|0.000000|0.000000|0.000000|1045
TCONS_00000010  XLOC_000007     AT1G01160|AT1G01160.2   =       q1:AT1G01160|AT1G01160.2|6|0.000000|0.000000|0.000000|1129
```   

the merged.tracking is the compact form of tmap file combined with the loci file. The last column has the gene, gene isoform, number of hits, FPKM, TPM, cov, and length all in that order.

Our last step is to create a count-table for the differential expression software "ballgown". A word of note about ballgown, ballgown requires that all of the files it will analyze be in their own parent directory, "ballgown", and furthermore each file is in its own subdirectory "ballgown/file_sub/file". Knowing we have four files, let's create the required structure for ballgown (run this):

```bash
mkdir ballgown
cd ballgown
mkdir athaliana_root_1 athaliana_root_2 athaliana_shoot_1 athaliana_shoot_2
```  

The reason for this will become obvious soon. Now we will use StringTie to make our count tables (do not run this!):  

```bash
module load stringtie/2.0.3
stringtie -e -B -p 16 athaliana_root_1.bam -G stringtie_merged.gtf -o ballgown/athaliana_root_1/athaliana_root_1.count

stringtie -e -B -p 16 athaliana_root_2.bam -G stringtie_merged.gtf -o ballgown/athaliana_root_2/athaliana_root_2.count

stringtie -e -B -p 16 athaliana_shoot_1.bam -G stringtie_merged.gtf -o ballgown/athaliana_shoot_1/athaliana_shoot_1.count

stringtie -e -B -p 16 athaliana_shoot_2.bam -G stringtie_merged.gtf -o ballgown/athaliana_shoot_2/athaliana_shoot_2.count
```  

Now we are ready to compile all of our code into a single script, which we call [transcript_assembly.sh](/04_assembly_n_quantify/transcript_assembly.sh)

```bash
module load gffread/0.9.12
module load stringtie/2.0.3

gffread TAIR10_GFF3_genes.gff -T -o athaliana_TAIR10_genes.gtf

#stringtie -p 16 ../03_align/athaliana_root_1.bam -G athaliana_TAIR10_genes.gtf -o athaliana_root_1.gtf
#stringtie -p 16 ../03_align/athaliana_root_2.bam -G athaliana_TAIR10_genes.gtf -o athaliana_root_2.gtf
#stringtie -p 16 ../03_align/athaliana_shoot_1.bam -G athaliana_TAIR10_genes.gtf -o athaliana_shoot_1.gtf
#stringtie -p 16 ../03_align/athaliana_shoot_2.bam -G athaliana_TAIR10_genes.gtf -o athaliana_shoot_2.gtf

stringtie --merge -p 16 -o stringtie_merged.gtf -G athaliana_TAIR10_genes.gtf mergelist.txt

module load gffcompare/0.10.4
gffcompare -r athaliana_TAIR10_genes.gtf -o merge stringtie_merged.gtf

mkdir ballgown
mkdir -p ballgown/athaliana_root_1 ballgown/athaliana_root_2 ballgown/athaliana_shoot_1 ballgown/athaliana_shoot_2


stringtie -e -B -p 16 ../03_align/athaliana_root_1.bam -G stringtie_merged.gtf -o ballgown/athaliana_root_1/athaliana_root_1.count
stringtie -e -B -p 16 ../03_align/athaliana_root_2.bam -G stringtie_merged.gtf -o ballgown/athaliana_root_2/athaliana_root_2.count
stringtie -e -B -p 16 ../03_align/athaliana_shoot_1.bam -G stringtie_merged.gtf -o ballgown/athaliana_shoot_1/athaliana_shoot_1.count
stringtie -e -B -p 16 ../03_align/athaliana_shoot_2.bam -G stringtie_merged.gtf -o ballgown/athaliana_shoot_2/athaliana_shoot_2.count
```

Once the stript has done executing you will see the following file structure in the **04_assembly_n_quantify** directory:
```
04_assembly_n_quantify/
├── athaliana_root_1.gtf
├── athaliana_root_2.gtf
├── athaliana_shoot_1.gtf
├── athaliana_shoot_2.gtf
├── athaliana_TAIR10_genes.gtf
├── ballgown
│   ├── athaliana_root_1
│   │   ├── athaliana_root_1.count
│   │   ├── e2t.ctab
│   │   ├── e_data.ctab
│   │   ├── i2t.ctab
│   │   ├── i_data.ctab
│   │   └── t_data.ctab
│   ├── athaliana_root_2
│   │   ├── athaliana_root_2.count
│   │   ├── e2t.ctab
│   │   ├── e_data.ctab
│   │   ├── i2t.ctab
│   │   ├── i_data.ctab
│   │   └── t_data.ctab
│   ├── athaliana_shoot_1
│   │   ├── athaliana_shoot_1.count
│   │   ├── e2t.ctab
│   │   ├── e_data.ctab
│   │   ├── i2t.ctab
│   │   ├── i_data.ctab
│   │   └── t_data.ctab
│   └── athaliana_shoot_2
│       ├── athaliana_shoot_2.count
│       ├── e2t.ctab
│       ├── e_data.ctab
│       ├── i2t.ctab
│       ├── i_data.ctab
│       └── t_data.ctab
├── merge.annotated.gtf
├── mergelist.txt
├── merge.loci
├── merge.stats
├── merge.stringtie_merged.gtf.refmap
├── merge.stringtie_merged.gtf.tmap
├── merge.tracking
├── stringtie_merged.gtf
├── TAIR10_GFF3_genes.gff
└── transcript_assembly.sh
```   


## 6. Differential expression analysis using ballgown  

For many organisms, many of the same genes are expressed in separate cell types, with a variety of phenotype differences a result of the specific isoforms a cell will use. Therefore, when performing a differential expression analysis from different parts of one organism (not one species, but a singular organism), it is wise to perform an isoform expression analysis alongside a standard differential expression analysis and combine the results (as we are doing here). We will only be performing the isoform expresion analysis. [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html) is a differential expression package for R via Bioconductor ideal for isoform expression analyses. Before beginning, you need to secure copy our **ballgown** directory from Xanadu to your local machine. 

To do this open a terminal window in your computer:
```bash
scp -r YOUR.USER.NAME@transfer.cam.uchc.edu:<PATH-TO-RNASeq_Model_Plant_Arabidopsis>/04_assembly_n_quantify/ballgown .
```  

   
Now we load [RStudio](https://www.rstudio.com/products/rstudio/download/) with administrator privileges (otherwise you cannot install packages!).

To begin we mush download and load the proper packages. Depending on the R version you are having in your local machine downloading instructions will differ. So please check the [biconductor website](https://www.bioconductor.org/) for appropriate instructions.   

Make sure you have the follwoing packages installed: 
```
devtools
RFLPtools
ballgown
genefilter
dplyr
ggplot2
gplots
RSkittleBrewer
```

Use Bioconductor version to install *ballgown, genefilter,dplyr,alyssafrazee/RSkittleBrewer* 
```
install.packages("RFLPtools")
install.packages("devtools")
install.packages("gplots")
#for R version 3.6 or later
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("alyssafrazee/RSkittleBrewer","ballgown", "genefilter", "dplyr", "ggplot2"))
```

Now we need to set our working directory to the directory which contains our "ballgown" folder. For me, this is:

```r
>dir <- "/Users/cbcuconn/RNASeq_Model_Plant_Arabidopsis/"

>setwd(dir)
>list.files()
```

You should see the "ballgown" folder after the list.files() command.

Let's have a look at the ballgown function:

```r
>help("ballgown")
```  

**constructor function for ballgown objects**
```
Usage  

ballgown(samples = NULL, dataDir = NULL, samplePattern = NULL,
  bamfiles = NULL, pData = NULL, verbose = TRUE, meas = "all")
Arguments

samples                 vector of file paths to folders containing sample-specific ballgown data (generated by tablemaker). If samples
                        is provided, dataDir and samplePattern are not used.
dataDir                 file path to top-level directory containing sample-specific folders with ballgown data in them. Only used if
                        samples is NULL.
samplePattern           regular expression identifying the subdirectories of\ dataDir containing data to be loaded into the ballgown
                        object (and only those subdirectories). Only used if samples is NULL.
bamfiles                optional vector of file paths to read alignment files for each sample. If provided, make sure to sort properly
                        (e.g., in the same order as samples). Default NULL.
pData                   optional data.frame with rows corresponding to samples and columns corresponding to phenotypic variables.
verbose                 if TRUE, print status messages and timing  Information as the object is constructed.
meas                    character vector containing either "all" or one or more of: "rcount", "ucount", "mrcount", "cov", "cov_sd",
                        "mcov", "mcov_sd", or "FPKM". The resulting ballgown object will only contain the specified expression
                        measurements, for the appropriate features. See vignette for which expression measurements are available for
                        which features. "all" creates the full object.
```  

Because of the structure of our ballgown directory, we may use dataDir = "ballgown", samplePattern = "athaliana", measure = "FPKM", and pData = some_type_of_phenotype_matrix.

We want all of the objects in our arguments to be in the same order as they are present in the ballgown directory. Therefore, we want our pData matrix to have two columns -- the first column being the samples as they appear in the ballgown directory, and the second being the phenotype of each sample in the column before it (root or shoot). Let's see the order of our sample files:  

```r
>list.files("ballgown/")

[1] "athaliana_root_1"  "athaliana_root_2"  "athaliana_shoot_1" "athaliana_shoot_2"
```  

Now we construct a 4x2 phenotype matrix with the first column being our samples in order and the second each sample's phenotype:










