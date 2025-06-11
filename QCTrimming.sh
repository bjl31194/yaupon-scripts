#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long
# Date : 3.24.2025
# Description: This script runs basic QC and trimming scripts on a set of fastq files.

## SLURM PARAMETERS ##
#SBATCH --job-name=fastqc_Ivo_full                            # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=32gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/log.%j			    # Standard output and error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/scratch/bjl31194/yaupon/wgs/plates234/trimmed_reads"
DATADIR="/scratch/bjl31194/yaupon/wgs/plates234"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#count reads
for file in $(ls *.fastq.gz); do echo $file &&
zcat $file | echo $((wc -l /4)); done > readcounts.txt

#load modules
ml FastQC/0.11.9-Java-11
ml Trim_Galore/0.6.7-GCCcore-11.2.0
ml cutadapt/4.9-GCCcore-12.3.0

# run FastQC
#fastqc /scratch/bjl31194/yaupon/wgs/plate1/raw_reads/*.fastq.gz --outdir $OUTDIR

# trim Illumina adapters
trim_galore --quality 20 --fastqc --illumina --retain_unpaired --paired $DATADIR/25055FL-01-01-01_S56_L006_R1_001.fastq.gz $DATADIR/25055FL-01-01-01_S56_L006_R2_001.fastq.gz