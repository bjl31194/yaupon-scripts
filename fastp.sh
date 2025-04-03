#!/bin/bash
#SBATCH --job-name=fastp_Ivo_array
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --time=24:00:00
#SBATCH --array=1-3
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

OUTDIR="/scratch/bjl31194/yaupon/wgs/plate1/trimmed_reads"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

# paths to reads
R1='/scratch/bjl31194/yaupon/wgs/plate1/raw_reads/'${name}'_L006_R1_001.fastq.gz'
R2='/scratch/bjl31194/yaupon/wgs/plate1/raw_reads/'${name}'_L006_R2_001.fastq.gz'

# load fastp
# https://github.com/OpenGene/fastp
ml fastp/0.23.2-GCC-11.3.0

# trim adapter sequences, remove reads shorter than 21nt, quality filtering, remove polyG tails
fastp -w 8 --dont_overwrite --in1 $R1 --in2 $R2 --out1 ${name}_R1_trimmed.fastq.gz --out2 ${name}_R2_trimmed.fastq.gz --unpaired1 ${name}_U.fastq.gz --unpaired2 ${name}_U.fastq.gz --failed_out ${name}_failed.fastq.gz -j ${name}.fastp_report.json -h ${name}.fastp_report.html
