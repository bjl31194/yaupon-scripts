#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long (adapted from Shreena Pradhan)
# Date : 06.01.2023
# Description: genotypes GVCF calls made by HaplotypeCaller
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=gatk_genotypegvcfs                                 # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=56gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/output/gatk_genotypegvcfs.%j.out			          # Standard output
#SBATCH --error=/scratch/bjl31194/output/gatk_genotypegvcfs.%j.err                # Error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


cd /scratch/bjl31194/yaupon

ml GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8

gatk --java-options "-Xmx64g" GenotypeGVCFs \
   -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta \
   -V gendb://yaupon_db \
   -O yaupon_cohort1_genotyped.vcf.gz
