#!/bin/bash
#SBATCH --job-name=snmf
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --time=15:00:00
#SBATCH --constraint=Genoa|Milan
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bjl31194@uga.edu

module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2

DATADIR="/scratch/bjl31194/yaupon/wgs/plates1234/vcf"

cd $DATADIR

R CMD BATCH /home/bjl31194/yaupon/yaupon-scripts/snmf.R