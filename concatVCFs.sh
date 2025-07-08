#!/bin/bash
#SBATCH --job-name=concatVCFs
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1234/vcf"

VCF_OUT="/scratch/bjl31194/yaupon/wgs/plates1234/vcf/Ivom_plates1234.vcf.gz"

# load modules
ml BCFtools/1.18-GCC-12.3.0

# move to the vcf directory
cd $DATADIR

# sort vcf files
bcftools sort Ilex_plates1234_1.vcf.gz -Ou -o Ilex_plates1234_1_sorted.bcf
bcftools sort Ilex_plates1234_2.vcf.gz -Ou -o Ilex_plates1234_2_sorted.bcf
bcftools sort Ilex_plates1234_3.vcf.gz -Ou -o Ilex_plates1234_3_sorted.bcf
bcftools sort Ilex_plates1234_4.vcf.gz -Ou -o Ilex_plates1234_4_sorted.bcf
bcftools sort Ilex_plates1234_5.vcf.gz -Ou -o Ilex_plates1234_5_sorted.bcf

# concatenate vcf files 
bcftools concat *_sorted.bcf -Oz --threads 8 -o Ilex_plates1234_merged.vcf.gz 

# optional filtering by sample id to remove decidua individuals (requires txt file with list of sample names to keep)
#bcftools view -Oz -S only_yaupon.txt Ivom_plate1_filter.vcf.gz > Ivom_plate1_sppfilter.vcf.gz