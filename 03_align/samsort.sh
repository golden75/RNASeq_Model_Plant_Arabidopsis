#!/bin/bash
#SBATCH --job-name=samsort
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


export TMPDIR=/labs/CBC/np_projects/RNASeq_Model_Plant_Arabidopsis/tmp

hostname

module load samtools/1.9

samtools sort -@ 6 -o athaliana_root_1.bam athaliana_root_1.sam
samtools sort -@ 6 -o athaliana_root_2.bam athaliana_root_2.sam
samtools sort -@ 6 -o athaliana_shoot_1.bam athaliana_shoot_1.sam
samtools sort -@ 6 -o athaliana_shoot_2.bam athaliana_shoot_2.sam


