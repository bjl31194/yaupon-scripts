#!/bin/bash
#SBATCH --job-name=bwa_ivom_array
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=80gb
#SBATCH --time=72:00:00
#SBATCH --array=14,23,29,68,74,83
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

OUTDIR="/scratch/bjl31194/yaupon/wgs/plate1/align"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

# name of assembly file
assembly='/scratch/bjl31194/yaupon/references/draft/I_vomitoria_GAF4_hap1_min50k.fa'

# paths to reads
R1='/scratch/bjl31194/yaupon/wgs/plate1/trimmed_reads/'${name}'_R1_trimmed.fastq.gz'
R2='/scratch/bjl31194/yaupon/wgs/plate1/trimmed_reads/'${name}'_R2_trimmed.fastq.gz'

# load samtools and bwa
ml SAMtools/1.16.1-GCC-11.3.0
ml BWA/0.7.17-GCCcore-11.3.0


​# map reads to new indexed references - hap1
bwa mem -t 32 $assembly $R1 $R2 | samtools view -@ 32 -O BAM | samtools sort -@ 32 -O BAM -o $name.Ivo.sorted.bam 

# index bam
samtools index $name.Ivo.sorted.bam
