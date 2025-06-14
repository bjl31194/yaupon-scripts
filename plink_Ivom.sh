#!/bin/bash
#SBATCH --job-name=plink_Ivom
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plate1/vcf"

VCF="/scratch/bjl31194/yaupon/wgs/plate1/vcf/Ivom_plate1_filter.vcf.gz"

# load modules
ml PLINK/2.0.0-a.6.9-gfbf-2023b
ml ADMIXTURE/1.3.0

# move to the vcf directory
cd $DATADIR

## Run plink to get .bed file and PCA ##

# identify prune sites
#plink --vcf $VCF --double-id --allow-extra-chr \
#--set-missing-var-ids @:# \
#--indep-pairwise 50 10 0.1 --out Ivom96

# linkage prune and create pca files
#plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
#--extract Ivom96.prune.in \
#--make-bed --pca --out Ivom96

## run ADMIXTURE ##

# generate input files
FILE=Ivom96
cd admixture

# Generate the input file in plink format
plink --vcf $VCF --make-bed --out $FILE --allow-extra-chr

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

# running ADMIXTURE for clusters size 2-5
for i in {2..5}
do
 admixture --cv $FILE.bed $i > log${i}.out
done

# yoink cross validation errors out of log files
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > $FILE.cv.error


