#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long (adapted from Shreena Pradhan)
# Date : 06.01.2023
# Description: creates merged genomicsdb from gvcf files
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=gatk_gdbi                                 # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/output/gatk_gdbi.%j.out			          # Standard output
#SBATCH --error=/scratch/bjl31194/output/gatk_gdbi.%j.err                # Error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


cd /scratch/bjl31194/yaupon/trimmed_reads

ml GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8

mkdir -p /lscratch/$SLURM_JOBID/tmp;

gatk --java-options "-Xmx8g -Xms8g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path /scratch/bjl31194/yaupon/yaupon_db \
       -L JYEU.hipmer.GA-F-4_scaffolds.intervals.list \
       --sample-name-map cohort1.sample_map \
       --tmp-dir /lscratch/$SLURM_JOBID/tmp \
       --reader-threads 4

cp -r /lscratch/${SLURM_JOB_ID} /scratch/bjl31194/output
