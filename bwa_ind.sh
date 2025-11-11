#!/bin/bash
#SBATCH --job-name=bwa_ind
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=7-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1234/align"

# name of assembly file (using haplotype 1)
assembly='/scratch/bjl31194/yaupon/references/draft/I_vomitoria_GAF4_hap1_min50k.fa'

# paths to reads
R1='/scratch/bjl31194/yaupon/wgs/plates1234/trimmed_reads/25055FL-03-01-63_S63_L007_R1_trimmed.fastq.gz'
R2='/scratch/bjl31194/yaupon/wgs/plates1234/trimmed_reads/25055FL-03-01-63_S63_L007_R2_trimmed.fastq.gz'

# load samtools and bwa
ml SAMtools/1.16.1-GCC-11.3.0
ml BWA/0.7.17-GCCcore-11.3.0


​# map reads to new indexed reference - hap1
bwa mem -t 32 $assembly $R1 $R2 | samtools view -@ 32 -O BAM | samtools sort -@ 32 -O BAM -o 25055FL-03-01-63_S63_L007.Ivo.sorted.bam 

# index bam
samtools index 25055FL-03-01-63_S63_L007.Ivo.sorted.bam