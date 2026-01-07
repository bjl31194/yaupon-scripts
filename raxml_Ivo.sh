#!/bin/bash
#SBATCH --job-name=raxml_Ilex
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=7-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# change sample names in vcf file
#ml BCFtools/1.18-GCC-12.3.0
#bcftools view -Ou Ivom_plate1_filter.vcf.gz | bcftools reheader -s names_plate1.txt -o Ivom_plate1_filter_names.vcf.gz

# set variables
OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/phylo"
VCF="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ilex_redrep_filter.vcf.gz"
PHYLIP="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ilex_redrep.min4.phy"

# load modules
ml RAxML-NG/1.2.2-GCC-13.2.0
ml Python/3.12.3-GCCcore-13.3.0

# move to the vcf directory
cd $OUTDIR

# build phylip matrix from vcf
python /home/bjl31194/yaupon/yaupon-scripts/vcf2phylip.py -i $VCF --output-folder $OUTDIR --output-prefix Ilex_redrep

# perform ML tree search and optimization
raxml-ng --all --bs-trees 1000 --msa $PHYLIP --model GTR+G
