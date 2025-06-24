#!/bin/bash
#SBATCH --job-name=raxml_Ivo
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

PHYLIP="/scratch/bjl31194/yaupon/wgs/plate1/vcf/Ivom_plate1_filter.min4.phy"

# load modules
ml RAxML-NG/1.2.0-GCC-12.3.0

# move to the vcf directory
cd $DATADIR

# perform ML tree search and optimization
raxml-ng --search1 --msa $PHYLIP --model GTR+G

