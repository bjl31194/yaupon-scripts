#!/bin/bash
#SBATCH --job-name=gemma_Ivom_wild
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --time=1-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

## set variables
OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/gwas"

SUBSET="Ivom_wild"

VCF="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ilex_plates1-5_names_filter_sexed_atlantic.vcf.gz"

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

## load modules
ml PLINK/2.0.0-a.6.20-gfbf-2024a
ml GEMMA/0.98.5-gfbf-2023b

## move to the proper directory 
cd $OUTDIR

## Make .bed file for inputting to GEMMA and filter data ##
# plink --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --nonfounders --set-missing-var-ids @:# \
# --maf 0.05 --geno 0.1 --mind 0.5 --snps-only \
# --make-bed --out Ilex_plates1-5_${SUBSET}

## attach phenotype data ##
plink --bfile ${SUBSET}_sexed --allow-no-sex --pheno ${SUBSET}_sex_phenotypes.txt \
--make-bed --out gemma_input_${SUBSET}

## make kinship matrix ##
gemma -bfile gemma_input_${SUBSET} -gk 1 -o RelMat_${SUBSET}

## run GEMMA (lmm=linear mixed model using kinship matrix, 2=likelihood ratio test) ##
gemma -bfile gemma_input_${SUBSET} -k output/RelMat_${SUBSET}.cXX.txt -lmm 2 -o GWAS_results_sex_${SUBSET}.lmm