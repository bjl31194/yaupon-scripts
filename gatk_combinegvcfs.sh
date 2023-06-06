#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long (adapted from Shreena Pradhan)
# Date : 06.06.2023
# Description: creates merged gvcf from per-sample gvcf files
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=gatk_combinegvcfs                               # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=64gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/bjl31194/output/gatk_combinegvcfs.%j.out			          # Standard output
#SBATCH --error=/scratch/bjl31194/output/gatk_combinegvcfs.%j.err                # Error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


cd /scratch/bjl31194/yaupon/trimmed_reads

ml GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8

gatk CombineGVCFs \
   -R /scratch/bjl31194/yaupon/references/JYEU.hipmer.GA-F-4_assembly.fasta \
   -V P2_A01.g.vcf.gz \
   -V P2_A03.g.vcf.gz \
   -V P2_A04.g.vcf.gz \
   -V P2_A05.g.vcf.gz \
   -V P2_B01.g.vcf.gz \
   -V P2_B02.g.vcf.gz \
   -V P2_B03.g.vcf.gz \
   -V P2_B04.g.vcf.gz \
   -V P2_B05.g.vcf.gz \
   -V P2_C01.g.vcf.gz \
   -V P2_C02.g.vcf.gz \
   -V P2_C03.g.vcf.gz \
   -V P2_C04.g.vcf.gz \
   -V P2_C05.g.vcf.gz \
   -V P2_D01.g.vcf.gz \
   -V P2_D02.g.vcf.gz \
   -V P2_D03.g.vcf.gz \
   -V P2_D04.g.vcf.gz \
   -V P2_D05.g.vcf.gz \
   -V P2_E01.g.vcf.gz \
   -V P2_E02.g.vcf.gz \
   -V P2_E03.g.vcf.gz \
   -V P2_E04.g.vcf.gz \
   -V P2_E05.g.vcf.gz \
   -V P2_F01.g.vcf.gz \
   -V P2_F02.g.vcf.gz \
   -V P2_F03.g.vcf.gz \
   -V P2_F04.g.vcf.gz \
   -V P2_F05.g.vcf.gz \
   -V P2_G01.g.vcf.gz \
   -V P2_G02.g.vcf.gz \
   -V P2_G03.g.vcf.gz \
   -V P2_G04.g.vcf.gz \
   -V P2_G05.g.vcf.gz \
   -V P2_H01.g.vcf.gz \
   -V P2_H02.g.vcf.gz \
   -V P2_H03.g.vcf.gz \
   -V P2_H04.g.vcf.gz \
   -V P2_H05.g.vcf.gz \
   -O yaupon_cohort1_combinedGVCF.g.vcf.gz
