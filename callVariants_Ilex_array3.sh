#!/bin/bash
#SBATCH --job-name=callVariants_Ilex_array20
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=7-00:00:00
#SBATCH --array=1-20
#SBATCH --mail-type=END,FAIL,ARRAY_TASKS
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

REGION=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/references/draft/filemap.txt)

# commands for making array filemap:
#cut -f1-2 I_vomitoria_GAF4_hap1_min50k.fa.fai > chrSize.txt
#awk '$1 = $1 FS "1"' chrSize.txt > chrSize_start.txt
#tr ' ' '\t' < chrSize_start.txt > chrSize_start_stop
#sort --random-sort chrSize_start_stop > contigs_shuffled
#split --numeric-suffixes=1 -n l/20 --additional-suffix='.txt' contigs_shuffled contigs

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

DATADIR="/scratch/bjl31194/yaupon/wgs/plates1234/align"

# name of assembly file
assembly='/scratch/bjl31194/yaupon/references/draft/I_vomitoria_GAF4_hap1_min50k.fa'

# load bcftools and samtools
ml BCFtools/1.18-GCC-12.3.0
ml SAMtools/1.16.1-GCC-11.3.0

# index reference
samtools faidx $assembly

# generate genotype likelihoods and call SNPs (no indels = -I option)
bcftools mpileup -d 100 -I -a AD,DP,SP,INFO/AD -Ou -f $assembly -b $DATADIR/bamlist.txt -R $REGION | bcftools call --threads 8 -mv -Oz -o $OUTDIR/Ilex_plates1-5_${SLURM_ARRAY_TASK_ID}.vcf.gz

