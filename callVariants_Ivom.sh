#!/bin/bash
#SBATCH --job-name=callVariants_Ivom
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=80gb
#SBATCH --time=72:00:00
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

# name of assembly file
assembly='/scratch/bjl31194/yaupon/references/draft/I_vomitoria_GAF4_hap1_min50k.fa'

# load bcftools and samtools
ml BCFtools/1.18-GCC-12.3.0
ml SAMtools/1.16.1-GCC-11.3.0

# index reference
samtools faidx $assembly

# pile up reads and call variants
bcftools mpileup -a AD,DP,SP -Oz -f $assembly /scratch/bjl31194/yaupon/wgs/plate1/align/*.sorted.bam | bcftools call -f GQ,GP \
-mO z -o $OUTDIR/Ivom_plate1.vcf.gz

