#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long
# Date : 06.18.2025
# Description: "So then I started blasting" -Frank
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=blast_TCS1                       # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out	                      # STDOUT
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error			                      # STDERR
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
SEQDIR="/scratch/bjl31194/yaupon/references/draft"
OUTDIR="/home/bjl31194/yaupon"
SUBJECT="/scratch/bjl31194/yaupon/references/draft/I_vomitoria_GAF4_hap1_min50k.fa"
QUERY="/scratch/bjl31194/yaupon/TCS1_CDS.fasta"

# download genomes
#curl http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz > ${SEQDIR}/Zm-B73-REFERENCE-NAM-5.0.fa.gz


# modules
ml BLAST+/2.16.0-gompi-2023b

#sed -n '1~4s/^@/>/p;2~4p' ${SEQDIR}/P2_A01.1.fq_trimmed.fq > ${SEQDIR}/P2_A01.1_trimmed.fasta

# make blast database
if [ ! -f ${SEQDIR}/I_vomitoria_GAF4_hap1_min50k_blastdb.ndb ]; then
	makeblastdb -dbtype nucl -in $SUBJECT -out ${SEQDIR}/I_vomitoria_GAF4_hap1_min50k_blastdb.ndb
fi

# blast sequence
blastn -num_threads 4 -query $QUERY -db ${SEQDIR}/I_vomitoria_GAF4_hap1_min50k_blastdb.ndb -out ${OUTDIR}/TCS1_blastresults.tsv -outfmt 6 -max_target_seqs 10
