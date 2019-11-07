#!/bin/bash
#SBATCH --job-name=sickle
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

module load sickle/1.33
 
sickle se -f ../01_raw_reads/SRR3498212.fastq -t sanger -o trimmed_SRR3498212.fastq -q 35 -l 45
sickle se -f ../01_raw_reads/SRR3498213.fastq -t sanger -o trimmed_SRR3498213.fastq -q 35 -l 45
sickle se -f ../01_raw_reads/SRR3498215.fastq -t sanger -o trimmed_SRR3498215.fastq -q 35 -l 45
sickle se -f ../01_raw_reads/SRR3498216.fastq -t sanger -o trimmed_SRR3498216.fastq -q 35 -l 45


