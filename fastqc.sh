#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long
# Date : 3.24.2025
# Description: This script runs basic QC and trimming scripts on a set of fastq files.

## SLURM PARAMETERS ##
#SBATCH --job-name=fastqc_Iv_plate1                            # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/log.%j			    # Standard output and error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/scratch/bjl31194/yaupon/wgs/plate1/trimmed_reads"
DATADIR="/work/gene8940/instructor_data"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#load modules
#ml FastQC/0.11.9-Java-11
ml Trimmomatic/0.39-Java-13

# run FastQC
#fastqc /scratch/bjl31194/yaupon/wgs/plate1/raw_reads/*.fastq.gz --outdir $OUTDIR

# trim Illumina adapters
java -jar trimmomatic/trimmomatic-0.39.jar \
    PE \
    -trimlog /scratch/bjl31194/yaupon/wgs/plate1/trimlog.txt \
    25055FL-01-01-01_S56_L006_R1_001.fastq.gz \
    25055FL-01-01-01_S56_L006_R2_001.fastq.gz \
    $OUTDIR/01_S56_R1_paired.fastq.gz \
    $OUTDIR/01_S56_R1_unpaired.fastq.gz \
    $OUTDIR/01_S56_R2_paired.fastq.gz \
    $OUTDIR/01_S56_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-SE.fa:2:30:10:2:True \
    LEADING:3 \
    TRAILING:3 \
    MINLEN:30
