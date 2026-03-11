#!/bin/bash
#SBATCH --job-name=alphafold
#SBATCH --partition=gpu_p
#SBATCH --gres=gpu:A100:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40gb
#SBATCH --time=7-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/alphafold"

cd $OUTDIR

ml AlphaFold/2.3.2-foss-2023a-CUDA-12.1.1

export ALPHAFOLD_DATA_DIR=/db/AlphaFold/2.3.2

alphafold --data_dir /db/AlphaFold/2.3.1 --output_dir . --fasta_paths ./AL-CW-6_ADH.fasta --max_template_date 2026-03-10
