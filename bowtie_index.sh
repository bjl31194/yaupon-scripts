#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long (adapted from Shreena Pradhan)
# Date : 05.31.2023
# Description: indexes a reference genome for aligning reads
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=indexing	                            # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/output/indexing.%j.out			          # Standard output
#SBATCH --error=/scratch/bjl31194/output/indexing.%j.err                # Error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/bjl31194/yaupon/references

#load modules
ml Bowtie2/2.4.1-GCC-8.3.0
ml picard/2.21.6-Java-11
ml SAMtools/1.10-iccifort-2019.5.281

#makes index files

bowtie2-build -f JYEU.hipmer.GA-F-4_assembly.fasta JYEU.hipmer.GA-F-4_assembly
samtools faidx JYEU.hipmer.GA-F-4_assembly.fasta
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=JYEU.hipmer.GA-F-4_assembly.fasta O=JYEU.hipmer.GA-F-4_assembly.dict
