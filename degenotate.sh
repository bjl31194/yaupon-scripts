#!/bin/bash
#SBATCH --job-name=mkado
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --time=3-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# set parameters
OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/mkado"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew"

REF="/scratch/bjl31194/yaupon/references/Ilex_vomitoria_var_GA_F_4_HAP1_V1_release/Ilex_vomitoria_var_GA_F_4/sequences/Ilex_vomitoria_var_GA_F_4.HAP1.mainGenome.fasta"
GFF="/scratch/bjl31194/yaupon/wgs/plates1-5/ann/Ilex_Hap1.filter.gff3"
VCF="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ivom1-5_names_newfilter.vcf.gz"
OUT="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/MC-IP-2.vcf.gz"
MAF=0.1
GENES="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/gene_IDs.txt"


cd $OUTDIR

# index vcf
#ml BCFtools/1.21-GCC-13.3.0
#bcftools index -t $VCF --threads 8 

conda init

# load conda environment
# conda activate degenotate
conda activate mkado

# run degenotate
# /home/bjl31194/.conda/envs/degenotate/bin/degenotate.py -a ${GFF} -g ${REF} -v ${VCF} -maf ${MAF} -e exclude.txt -u outgroup_Ipa.txt -o ${OUTDIR} -sfs --overwrite

# run mkado
/home/bjl31194/.conda/envs/mkado/bin/mkado vcf --vcf ${VCF} --outgroup-vcf ${OUT} --ref ${REF} --gff ${GFF} --per-gene --workers 8 --format tsv --verbose -a 
