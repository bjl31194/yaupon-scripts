#!/bin/bash
#SBATCH --job-name=calc_pairwiseFst
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=7-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1234/vcf"
SCRIPTS="/home/bjl31194/yaupon/yaupon-scripts"
ml R/4.4.2-gfbf-2024a

cd $DATADIR

R CMD INSTALL --library=/home/bjl31194/Rlibs adegenet_2.1.11.tar.gz
R CMD INSTALL --library=/home/bjl31194/Rlibs hierfstat_0.5-11.tar.gz

R CMD BATCH $SCRIPTS/pairwiseFst.R