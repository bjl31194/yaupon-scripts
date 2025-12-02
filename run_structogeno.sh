#!/bin/bash
#SBATCH --job-name=run_structogeno
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bjl31194@uga.edu

DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcf"

cd $DATADIR

module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2 # this will load R package LEA

R CMD BATCH /home/bjl31194/yaupon/yaupon-scripts/structogeno.R