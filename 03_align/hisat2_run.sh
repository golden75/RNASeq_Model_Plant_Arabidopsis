#!/bin/bash
#SBATCH --job-name=hisat2_run
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=50G
#SBATCH --partition=himem2
#SBATCH --qos=himem
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


export TMPDIR=/labs/CBC/np_projects/RNASeq_Model_Plant_Arabidopsis/tmp

hostname

module load hisat2/2.1.0

hisat2 -p 6 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q ../02_trimmed_reads/trimmed_SRR3498212.fastq  -S athaliana_root_1.sam 
hisat2 -p 6 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q ../02_trimmed_reads/trimmed_SRR3498213.fastq -S athaliana_root_2.sam
hisat2 -p 6 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q ../02_trimmed_reads/trimmed_SRR3498215.fastq -S athaliana_shoot_1.sam
hisat2 -p 6 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/Athaliana_HISAT2/thaliana -q ../02_trimmed_reads/trimmed_SRR3498216.fastq -S athaliana_shoot_2.sam




