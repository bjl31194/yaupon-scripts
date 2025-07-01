#!/bin/bash
#SBATCH --job-name=filterVariants_Ivom
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plate1/vcf"

VCF_IN="/scratch/bjl31194/yaupon/wgs/plate1/vcf/Ivom_plate1.vcf.gz"
VCF_OUT="/scratch/bjl31194/yaupon/wgs/plate1/vcf/Ivom_plate1_filter.vcf.gz"

# set filters
MAF=0.01 # median was 0.02, 1st quantile was 0.01. keeping this filter weak for now
MISS=0.9 # mean was 0.95 (only missing 5%)
QUAL=30 # 0.01% error rate or lower
MIN_DEPTH=8 # 1st quantile 7.63, mean was 9
MAX_DEPTH=25 # went with mean depth x 2 + a little extra judging from where most variants fell on histogram

# load modules
ml VCFtools/0.1.16-GCC-13.3.0
#ml BCFtools/1.18-GCC-12.3.0

# move to the vcf directory
cd $DATADIR

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT

# optional filtering by sample id to remove decidua individuals (requires txt file with list of sample names to keep)
#bcftools view -Oz -S only_yaupon.txt Ivom_plate1_filter.vcf.gz > Ivom_plate1_sppfilter.vcf.gz