#!/bin/bash
#SBATCH --job-name=quality_check
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


export TMPDIR=/labs/CBC/np_projects/RNASeq_Model_Plant_Arabidopsis/tmp

hostname


module load fastqc/0.11.7
module load MultiQC/1.7

fastqc trimmed_SRR3498212.fastq
fastqc trimmed_SRR3498213.fastq
fastqc trimmed_SRR3498215.fastq
fastqc trimmed_SRR3498216.fastq

multiqc -f -n trimmed trimmed*


