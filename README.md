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


The sort function converts SAM files to BAM automatically. Therefore, we can cut through most of these options and do a simple "samtools sort -o <output.bam> <inupt.sam>. Let's write our script:





