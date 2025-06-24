#!/bin/bash
#SBATCH --job-name=raxml_Ivo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# change sample names in vcf file
#ml BCFtools/1.18-GCC-12.3.0
#bcftools view -Ou Ivom_plate1_filter.vcf.gz | bcftools reheader -s names_plate1.txt -o Ivom_plate1_filter_names.vcf.gz

# set variables
DATADIR="/scratch/bjl31194/yaupon/wgs/plate1/vcf"
VCF="/scratch/bjl31194/yaupon/wgs/plate1/vcf/Ivom_plate1_filter_names.vcf.gz"
PHYLIP="/scratch/bjl31194/yaupon/wgs/plate1/vcf/Ivom_plate1_filter_names.min4.phy"

# load modules
ml RAxML-NG/1.2.0-GCC-12.3.0
ml Python/3.12.3-GCCcore-13.3.0

# move to the vcf directory
cd $DATADIR

# build phylip matrix from vcf
python /home/bjl31194/yaupon/yaupon-scripts/vcf2phylip.py -i $VCF --output-folder $DATADIR

# perform ML tree search and optimization
raxml-ng --search1 --msa $PHYLIP --model GTR+G