#!/bin/bash
#SBATCH --job-name=callVariants_Ivom
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

OUTDIR="/scratch/bjl31194/yaupon/wgs/plate1/vcf"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

DATADIR="/scratch/bjl31194/yaupon/wgs/plate1/align"

# name of assembly file
assembly='/scratch/bjl31194/yaupon/references/draft/I_vomitoria_GAF4_hap1_min50k.fa'

# load bcftools and samtools
ml BCFtools/1.18-GCC-12.3.0
ml SAMtools/1.16.1-GCC-11.3.0

# index reference
samtools faidx $assembly

# generate genotype likelihoods and call SNPs (no indels = -I option)
bcftools mpileup -I -a AD,DP,SP -Ou -f $assembly $DATADIR/*.sorted.bam | bcftools call --threads 32 -mv -Oz -o $OUTDIR/Ivom_plate1.vcf.gz

# test with 3 samples - looks good, took 2 hours to run tho
#bcftools mpileup -a AD,DP,SP -Ou -f $assembly $DATADIR/25055FL-01-01-01_S56.Ivo.sorted.bam $DATADIR/25055FL-01-01-02_S57.Ivo.sorted.bam $DATADIR/25055FL-01-01-03_S58.Ivo.sorted.bam | bcftools call -f GQ,GP -mO z -o $OUTDIR/Ivom_plate1.vcf.gz