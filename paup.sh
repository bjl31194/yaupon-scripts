#!/bin/bash
#SBATCH --job-name=paup_svdq
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=7-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/paup"

cd $OUTDIR

module load PAUP/4.0a168-centos64

paup Ilex_spp.nex -L Ilex_spp.log -n

svdquartets evalquartets=all nthreads=8 bootstrap;