#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long
# Date : ..2021
# Description: :)
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=JOB	                                # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=1:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/bjl31194/log.%j			    # Standard output and error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/work/gene8940/bjl31194/"
DATADIR="/work/gene8940/instructor_data"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#load modules
