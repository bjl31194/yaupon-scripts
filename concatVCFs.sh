#!/bin/bash
#SBATCH --job-name=concatVCFs
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew"

VCF_OUT="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ilex_plates1-5.vcf.gz"

# load modules
ml BCFtools/1.18-GCC-12.3.0

# move to the vcf directory
cd $DATADIR

# sort vcf files
bcftools sort Ilex_plates1-5_1.vcf.gz -Ou -o Ilex_plates1-5_1_sorted.bcf
bcftools sort Ilex_plates1-5_2.vcf.gz -Ou -o Ilex_plates1-5_2_sorted.bcf
bcftools sort Ilex_plates1-5_3.vcf.gz -Ou -o Ilex_plates1-5_3_sorted.bcf
bcftools sort Ilex_plates1-5_4.vcf.gz -Ou -o Ilex_plates1-5_4_sorted.bcf
bcftools sort Ilex_plates1-5_5.vcf.gz -Ou -o Ilex_plates1-5_5_sorted.bcf
bcftools sort Ilex_plates1-5_6.vcf.gz -Ou -o Ilex_plates1-5_6_sorted.bcf
bcftools sort Ilex_plates1-5_7.vcf.gz -Ou -o Ilex_plates1-5_7_sorted.bcf
bcftools sort Ilex_plates1-5_8.vcf.gz -Ou -o Ilex_plates1-5_8_sorted.bcf
bcftools sort Ilex_plates1-5_9.vcf.gz -Ou -o Ilex_plates1-5_9_sorted.bcf
bcftools sort Ilex_plates1-5_10.vcf.gz -Ou -o Ilex_plates1-5_10_sorted.bcf
bcftools sort Ilex_plates1-5_11.vcf.gz -Ou -o Ilex_plates1-5_11_sorted.bcf
bcftools sort Ilex_plates1-5_12.vcf.gz -Ou -o Ilex_plates1-5_12_sorted.bcf
bcftools sort Ilex_plates1-5_13.vcf.gz -Ou -o Ilex_plates1-5_13_sorted.bcf
bcftools sort Ilex_plates1-5_14.vcf.gz -Ou -o Ilex_plates1-5_14_sorted.bcf
bcftools sort Ilex_plates1-5_15.vcf.gz -Ou -o Ilex_plates1-5_15_sorted.bcf
bcftools sort Ilex_plates1-5_16.vcf.gz -Ou -o Ilex_plates1-5_16_sorted.bcf
bcftools sort Ilex_plates1-5_17.vcf.gz -Ou -o Ilex_plates1-5_17_sorted.bcf
bcftools sort Ilex_plates1-5_18.vcf.gz -Ou -o Ilex_plates1-5_18_sorted.bcf
bcftools sort Ilex_plates1-5_19.vcf.gz -Ou -o Ilex_plates1-5_19_sorted.bcf
bcftools sort Ilex_plates1-5_20.vcf.gz -Ou -o Ilex_plates1-5_20_sorted.bcf


# concatenate vcf files 
bcftools concat *_sorted.bcf -Oz --threads 8 -o Ilex_plates1-5_merged.vcf.gz 

# optional filtering by sample id to remove decidua individuals (requires txt file with list of sample names to keep)
#bcftools view -Oz -S only_yaupon.txt Ivom_plate1_filter.vcf.gz > Ivom_plate1_sppfilter.vcf.gz