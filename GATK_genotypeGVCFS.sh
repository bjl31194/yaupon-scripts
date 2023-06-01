#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long (adapted from Shreena Pradhan)
# Date : 06.01.2023
# Description: genotypes GVCF calls made by hapcaller
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=gatk_hapcaller                                 # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/output/gatk_genotypegvcfs.%j.out			          # Standard output
#SBATCH --error=/scratch/bjl31194/output/gatk_genotypegvcfs.%j.err                # Error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


cd /scratch/bjl31194/yaupon/trimmed_reads

ml GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8

gatk --java-options "-Xmx10g" GenotypeGVCFs \
   -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta \
   -V P2_A01.g.vcf.gz \
   -O P2_A01_genotyped.vcf.gz
