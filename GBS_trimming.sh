#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long (adapted from Shreena Pradhan)
# Date : 4.28.2021
# Description: trims GBS reads using fastx_trimmer
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=trimming	                            # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=trimming.%j.out			          # Standard output
#SBATCH --error=trimming.%j.err                # Error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/scratch/bjl31194/yaupon/sequences/Fastqc_results/trimmed_reads"
DATADIR="/home/bjl31194/yaupon/sequences/Yaupon_barcode_filter_reads"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#load modules
module load FASTX-Toolkit/0.0.14-GCCcore-8.3.0

#trimming

cd ${OUTDIR}

for j in ${DATADIR}/*.1.fq;

do

fastx_trimmer -Q 33 -f 6 -l 86 -i $j -o ${j}_trimmed.fq

done

for k in ${DATADIR}/*.2.fq;

do
fastx_trimmer -Q 33 -f 6 -l 96 -i $k -o ${k}_trimmed.fq

done
