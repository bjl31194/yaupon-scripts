#!/bin/bash
#SBATCH --job-name=index_Ivom
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

#path to reference
REF='/scratch/bjl31194/yaupon/references/draft/ '

#generate index files in same directory as reference
bwa index $REF 