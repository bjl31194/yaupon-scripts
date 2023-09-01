#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long (adapted from Shreena Pradhan)
# Date : 08.31.2023
# Description: filters variants using VariantFiltration
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=gatk_variantfiltration                                 # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=56gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/output/gatk_variantfiltration.%j.out			          # Standard output
#SBATCH --error=/scratch/bjl31194/output/gatk_variantfiltration.%j.err                # Error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


cd /scratch/bjl31194/yaupon/vcf

ml GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8

gatk SelectVariants \
   -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta \
   -V yaupon_cohort1_genotyped.vcf \
   --select-type-to-include SNP \
   --restrict-alleles-to BIALLELIC \
   -O cohort1_biallelic_SNPs.vcf \

gatk VariantFiltration \
    -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta \
    -V cohort1_biallelic_SNPs.vcf \
    -O cohort1_biallelic_QD08_SNPs.vcf \
    --filter-name "QD08" \
    --filter-expression "QD > 8.00" \

gatk VariantsToTable \
     -V cohort1_QD10_SNPs.vcf.gz \
     -F CHROM -F POS -F TYPE -F QD -F MQ -GF AD \
     -O yaupon_cohort1_genotyped_filtered.table
