#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long
# Date : 3.18.2025
# Description: This script runs fastQC on a set of fastq files.
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=fastqc_Iv_plate1                            # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=1:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/log.%j			    # Standard output and error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/scratch/bjl31194/yaupon/wgs/plate1/fastqc"
DATADIR="/work/gene8940/instructor_data"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#load modules
ml FastQC/0.11.9-Java-11

# run FastQC
fastqc /scratch/bjl31194/yaupon/wgs/plate1/raw_reads/*.fastq.gz --outdir $OUTDIR