#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long (adapted from Shreena Pradhan)
# Date : 05.31.2023
# Description: aligns reads to a reference using bowtie
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=indexing	                            # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/output/bowtie.%j.out			          # Standard output
#SBATCH --error=/scratch/bjl31194/output/bowtie.%j.err                # Error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


cd /scratch/bjl31194/yaupon/trimmed_reads

ml Bowtie2/2.4.1-GCC-8.3.0
ml picard/2.21.6-Java-11
ml SAMtools/1.10-iccifort-2019.5.281

bowtie2 -x JYEU.hipmer.GA-F-4_assembly -1 pqX.1.fq.fq -2 pqX.2.fq.fq -p 4 -S pqX.sam
samtools view -@ 4 -bS pqX.sam > pqX.bam
java -jar  $EBROOTPICARD/picard.jar ValidateSamFile I=pqX.bam
java -jar  $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=pqX.bam O=pqX.Gr.bam RGLB=Whatever RGPU=Whatever RGPL=illumina RGSM=pqX
samtools sort pqX.Gr.bam -o pqX.Gr.sorted.bam
samtools index pqX.Gr.sorted.bam
