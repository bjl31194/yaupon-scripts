#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long (adapted from Shreena Pradhan)
# Date : 05.31.2023
# Description: calls variants using UnifiedGenotyper in GATK
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=gatk_ug                                 # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=64gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/output/gatk.%j.out			          # Standard output
#SBATCH --error=/scratch/bjl31194/output/gatk.%j.err                # Error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


cd /scratch/bjl31194/yaupon/trimmed_reads

ml GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta \
-I P2_A01.Gr.sorted.bam \
-I P2_A03.Gr.sorted.bam \
-I P2_A05.Gr.sorted.bam \
-I P2_B01.Gr.sorted.bam \
-I P2_B02.Gr.sorted.bam \
-I P2_B03.Gr.sorted.bam \
-I P2_B04.Gr.sorted.bam \
-I P2_B05.Gr.sorted.bam \
-I P2_C01.Gr.sorted.bam \
-I P2_C02.Gr.sorted.bam \
-I P2_C03.Gr.sorted.bam \
-I P2_C04.Gr.sorted.bam \
-I P2_C05.Gr.sorted.bam \
-I P2_D01.Gr.sorted.bam \
-I P2_D02.Gr.sorted.bam \
-I P2_D03.Gr.sorted.bam \
-I P2_D04.Gr.sorted.bam \
-I P2_D05.Gr.sorted.bam \
-I P2_E01.Gr.sorted.bam \
-I P2_E02.Gr.sorted.bam \
-I P2_E03.Gr.sorted.bam \
-I P2_E04.Gr.sorted.bam \
-I P2_E05.Gr.sorted.bam \
-I P2_F01.Gr.sorted.bam \
-I P2_F02.Gr.sorted.bam \
-I P2_F03.Gr.sorted.bam \
-I P2_F04.Gr.sorted.bam \
-I P2_F05.Gr.sorted.bam \
-I P2_G01.Gr.sorted.bam \
-I P2_G02.Gr.sorted.bam \
-I P2_G03.Gr.sorted.bam \
-I P2_G04.Gr.sorted.bam \
-I P2_G05.Gr.sorted.bam \
-I P2_H01.Gr.sorted.bam \
-I P2_H02.Gr.sorted.bam \
-I P2_H03.Gr.sorted.bam \
-I P2_H04.Gr.sorted.bam \
-I P2_H05.Gr.sorted.bam \
-dcov 1000 -glm BOTH -o yaupon_cohort1_raw_SNPs.vcf
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta -o yaupon_cohort1_biallelic_raw_SNPs.vcf --variant yaupon_cohort1_raw_SNPs.vcf -restrictAllelesTo BIALLELIC
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta -V yaupon_cohort1_biallelic_raw_SNPs.vcf --filterExpression "QD <= 10.00" --filterName "QD_10" -o yaupon_cohort1_QD10_biallelic_raw_SNPs.vcf
