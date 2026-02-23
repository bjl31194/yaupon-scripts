#!/bin/bash
#SBATCH --job-name=orthofinder_Hannuus
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=3-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/orthofinder"

cd $DATADIR

# load modules
ml OrthoFinder/3.1.0-foss-2023a

# run orthofinder (very hard)
orthofinder -f .
