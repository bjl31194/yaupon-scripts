#!/bin/bash
#SBATCH --job-name=raxml_Ivo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# change sample names in vcf file
#ml BCFtools/1.18-GCC-12.3.0
#bcftools view -Ou Ivom_plate1_filter.vcf.gz | bcftools reheader -s names_plate1.txt -o Ivom_plate1_filter_names.vcf.gz

# set variables
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1234/vcf"
VCF="/scratch/bjl31194/yaupon/wgs/plates1234/vcf/Ilex_plates1234_filtered.vcf.gz"
PHYLIP="/scratch/bjl31194/yaupon/wgs/plate1/vcf/Ilex_plates1234_filtered.min4.phy"

# load modules
ml RAxML-NG/1.2.0-GCC-12.3.0
#ml Python/3.12.3-GCCcore-13.3.0

# move to the vcf directory
cd $DATADIR

# build phylip matrix from vcf
#python /home/bjl31194/yaupon/yaupon-scripts/vcf2phylip.py -i $VCF --output-folder $DATADIR

# perform ML tree search and optimization
raxml-ng --all --bs-trees 100 --msa $PHYLIP --model GTR+G
