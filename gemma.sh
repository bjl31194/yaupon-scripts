#!/bin/bash
#SBATCH --job-name=gemma_Ilex
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --time=1-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/gwas"

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/gwas"

SUBSET="gulf"

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load modules
ml PLINK/2.0.0-a.6.20-gfbf-2024a
ml GEMMA/0.98.5-gfbf-2023b

## move to the proper directory
cd $DATADIR

## make kinship matrix
gemma -bfile gemma_input_${SUBSET} -gk 1 -o RelMat_${SUBSET}

## run GEMMA (lmm=linear mixed model using kinship matrix, 2=likelihood ratio test)
gemma -bfile gemma_input_${SUBSET} -k output/RelMat_${SUBSET}.cXX.txt -lmm 2 -o GWAS_results_sex_${SUBSET}.lmm