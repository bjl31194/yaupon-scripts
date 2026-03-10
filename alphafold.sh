#!/bin/bash
#SBATCH --job-name=alphafold
#SBATCH --partition=batch
#SBATCH --partition=gpu_p
#SBATCH --gres=gpu:L4:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=700gb
#SBATCH --time=7-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/alphafold"

cd $OUTDIR

ml CUDA/12.1.1
ml AlphaFold/2.3.2-foss-2023a-CUDA-12.1.1

alphafold --fasta_paths ./AL-CW-6_ADH.fasta --output_dir . --max_template_date 3000-01-01