#!/bin/bash
#SBATCH --job-name=mapping_stats
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=1-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/align_hap2"

# load samtools
ml SAMtools/1.16.1-GCC-11.3.0

cd $OUTDIR 

( for i in $OUTDIR/*.bam ; do samtools flagstat $i | awk -F "[(|%]" 'NR == 7 {print $2}' ; done) > mapping_percentages.txt