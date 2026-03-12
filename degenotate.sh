#!/bin/bash
#SBATCH --job-name=degenotate
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=3-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# set parameters
OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/degenotate"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew"

REF="/scratch/bjl31194/yaupon/references/Ilex_vomitoria_var_GA_F_4_HAP1_V1_release/Ilex_vomitoria_var_GA_F_4/sequences/Ilex_vomitoria_var_GA_F_4.HAP1.mainGenome.fasta"
GFF="/scratch/bjl31194/yaupon/wgs/plates1-5/ann/Ilex_Hap1.filter.gff3"
VCF="Ilex1-5_names.vcf.gz"
MAF=0.1


cd $DATADIR

conda init

# load conda environment
conda activate degenotate

# run degenotate
degenotate.py -a ${GFF} -g ${REF} -v ${VCF} -maf ${MAF} -e exclude.txt -u outgroup_Ipa.txt -o ${OUTDIR} -sfs
