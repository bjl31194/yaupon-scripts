#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long
# Date : 05.24.2023
# Description: "So then I started blasting" -Frank
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=blast_homologs                       # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=8		                            # Number of cores per task
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=~/yaupon/blast.%j.out			                      # STDOUT
#SBATCH --error=~/yaupon/blast.%j.err			                      # STDERR
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
SEQDIR="/home/bjl31194/yaupon/sequences/Yaupon_barcode_filter_reads"
OUTDIR="/home/bjl31194/yaupon"

#download genomes

#curl http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz > ${SEQDIR}/Zm-B73-REFERENCE-NAM-5.0.fa.gz


#modules
ml BLAST+/2.11.0-gompi-2019b
module list
​
sed -n '1~4s/^@/>/p;2~4p' ${SEQDIR}/P2_A01.1.fq_trimmed.fq > ${SEQDIR}/P2_A01.1_trimmed.fasta

#make blast databases
if [ ! -f ${SEQDIR}/P2_A01.1_blastdb.ndb ]; then
	makeblastdb -dbtype nucl -in ${SEQDIR}/P2_A01.1_trimmed.fasta -out ${SEQDIR}/P2_A01.1_blastdb.ndb
fi

#blast sequences to maize db
blastn -num_threads 8 -query /home/bjl31194/yaupon/gene.fna -db ${SEQDIR}/P2_A01.1_blastdb.ndb -out ${OUTDIR}/TCS1_blastresults.tsv -outfmt 6 -max_target_seqs 10
