#!/bin/bash
#SBATCH --job-name=transcript_assembly.sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


export TMPDIR=/labs/CBC/np_projects/RNASeq_Model_Plant_Arabidopsis/tmp

hostname
module load gffread/0.9.12
module load stringtie/2.0.3

#gffread TAIR10_GFF3_genes.gff -T -o athaliana_TAIR10_genes.gtf

#stringtie -p 16 ../03_align/athaliana_root_1.bam -G athaliana_TAIR10_genes.gtf -o athaliana_root_1.gtf
#stringtie -p 16 ../03_align/athaliana_root_2.bam -G athaliana_TAIR10_genes.gtf -o athaliana_root_2.gtf
#stringtie -p 16 ../03_align/athaliana_shoot_1.bam -G athaliana_TAIR10_genes.gtf -o athaliana_shoot_1.gtf
#stringtie -p 16 ../03_align/athaliana_shoot_2.bam -G athaliana_TAIR10_genes.gtf -o athaliana_shoot_2.gtf

#stringtie --merge -p 16 -o stringtie_merged.gtf -G athaliana_TAIR10_genes.gtf mergelist.txt 
#module load gffcompare/0.10.4
#gffcompare -r athaliana_TAIR10_genes.gtf -o merge stringtie_merged.gtf

mkdir ballgown
mkdir -p ballgown/athaliana_root_1 ballgown/athaliana_root_2 ballgown/athaliana_shoot_1 ballgown/athaliana_shoot_2


stringtie -e -B -p 16 ../03_align/athaliana_root_1.bam -G stringtie_merged.gtf -o ballgown/athaliana_root_1/athaliana_root_1.count
stringtie -e -B -p 16 ../03_align/athaliana_root_2.bam -G stringtie_merged.gtf -o ballgown/athaliana_root_2/athaliana_root_2.count
stringtie -e -B -p 16 ../03_align/athaliana_shoot_1.bam -G stringtie_merged.gtf -o ballgown/athaliana_shoot_1/athaliana_shoot_1.count
stringtie -e -B -p 16 ../03_align/athaliana_shoot_2.bam -G stringtie_merged.gtf -o ballgown/athaliana_shoot_2/athaliana_shoot_2.count




