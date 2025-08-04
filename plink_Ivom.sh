#!/bin/bash
#SBATCH --job-name=strct_th_Ilex
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=7-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1234/vcf/structure"

VCF="/scratch/bjl31194/yaupon/wgs/plates1234/vcf/Ilex_plates1234_filtered.vcf.gz"

STRUCT_IN="/scratch/bjl31194/yaupon/wgs/plates1234/vcf/structure/Ilex1234forStructure.recode.strct_in"

# load modules
#ml PLINK/2.0.0-a.6.9-gfbf-2023b
#ml ADMIXTURE/1.3.0
ml Structure/2.3.4-GCC-12.3.0
ml Structure_threader/1.3.10-foss-2023a

# move to the proper directory
cd $DATADIR

## Run plink to get .bed file and PCA ##

# identify prune sites
#plink --vcf $VCF --double-id --allow-extra-chr \
#--set-missing-var-ids @:# \
#--indep-pairwise 50 10 0.1 --out Ivom384

# linkage prune and create pca files
#plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
#--extract Ivom384.prune.in \
#--make-bed --pca --out Ivom384

# generate structure input file
#plink --bfile Ivom384 --allow-extra-chr --recode structure --out Ivom384forStructure

## run ADMIXTURE ##

# generate input files
#FILE=Ilex384

# Generate the input file in plink format
#plink --vcf $VCF --double-id --make-bed --out $FILE --allow-extra-chr

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
#awk '{$1="0";print $0}' $FILE.bim > $FILE.bim.tmp
#mv $FILE.bim.tmp $FILE.bim

# running ADMIXTURE for clusters size 3-8
#for i in {3..8}
#do
#    admixture --cv $FILE.bed $i > log${i}.out
#done

# yoink cross validation errors out of log files
#awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > $FILE.cv.error

## STRUCTURE - for running on cluster ##

structure_threader run -K 5 -R 4 -i $STRUCT_IN -o $DATADIR -t 32 --ind indfile.csv -st /apps/eb/Structure/2.3.4-GCC-12.3.0/bin/structure

## other misc scripts ##

# run pong locally for ADMIXTURE visualization
# use Q matrix files from ADMIXTURE output
#pong -m pong_filemap.txt -i ind2pop.txt
