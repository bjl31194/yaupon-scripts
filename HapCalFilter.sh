#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long (adapted from Shreena Pradhan)
# Date : 08.31.2023
# Description: filters variants using HC
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=hapcalfilter                                # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=8		                            # Number of cores per task
#SBATCH --mem=32gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/output/gatk_hapfil.%j.out			          # Standard output
#SBATCH --error=/scratch/bjl31194/output/gatk_hapfil.%j.err                # Error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/bjl31194/yaupon/vcf

ml GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants \ -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta -o outputhapF.vcf --variant yaupon_cohort1_genotyped.vcf -restrictAllelesTo BIALLELIC
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta -V outputhapF.vcf --filterExpression "QD <= 10.00" --filterName "QD_10" -o cohort1_biallelicQD10_fromHC_SNPs.vcf
