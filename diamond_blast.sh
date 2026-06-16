#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long
# Date : 06.18.2025
# Description: "So then I started blasting" -Frank
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=diamond                                  # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=32		                            # Number of cores per task
#SBATCH --mem=32gb			                                # Total memory for job
#SBATCH --time=7-00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out	                      # STDOUT
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error			                      # STDERR
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

# modules
ml BLAST+/2.17.0-gompi-2025a
ml ncbiblastdb/20260501
ml DIAMOND/2.1.23-GCC-13.3.0

# Variables
OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/diamond"
prot_file=/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/diamond/Ilex_Hap1_peptides.filter.fasta
out_file=Ivo_Hap1_protein_ann_uniprot
db=/db/ncbiblast/20260501/swissprot

cd $OUTDIR

# Run DIAMOND
diamond blastp --threads $SLURM_CPUS_PER_TASK --max-target-seqs 5 --evalue 0.001 \
               --db $db --query ${prot_file} --outfmt 6 \
               --out gene_annotations/${out_file}.uniprot.txt