#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long (adapted from Shreena Pradhan)
# Date : 05.31.2023
# Description: calls variants using UnifiedGenotyper in GATK
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=gatk                                 # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/output/gatk.%j.out			          # Standard output
#SBATCH --error=/scratch/bjl31194/output/gatk.%j.err                # Error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


cd /scratch/bjl31194/yaupon/trimmed_reads

ml GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta -I P2_A01.Gr.sorted.bam -I P2_A03.Gr.sorted.bam -dcov 1000 -glm BOTH -o raw_SNPs.vcf
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta -o biallelic_raw_SNPs.vcf --variant raw_SNPs.vcf -restrictAllelesTo BIALLELIC
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta -V biallelic_raw_SNPs.vcf --filterExpression "QD <= 10.00" --filterName "QD_10" -o QD10_biallelic_raw_SNPs.vcf
