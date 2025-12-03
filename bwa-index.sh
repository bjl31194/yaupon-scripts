#!/bin/bash
#SBATCH --job-name=index_Ivom_hap2
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
REF1='/scratch/bjl31194/yaupon/references/Ilex_vomitoria_var_GA_F_4_HAP1_V1_release/Ilex_vomitoria_var_GA_F_4/sequences/Ilex_vomitoria_var_GA_F_4.HAP1.mainGenome.fasta'
REF2='/scratch/bjl31194/yaupon/references/Ilex_vomitoria_var_GA_F_4_HAP2_V1_release/Ilex_vomitoria_var_GA_F_4/sequences/Ilex_vomitoria_var_GA_F_4.HAP2.mainGenome.fasta'
# load module
ml BWA/0.7.17-GCCcore-11.3.0

#generate index files in same directory as reference
# bwa index $REF1

bwa index $REF2