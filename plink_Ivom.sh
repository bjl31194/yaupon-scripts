#!/bin/bash
#SBATCH --job-name=plink_Ilex
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=1-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/gwas"

VCF="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/gwas/Ilex_plates1-5_names_filter_texas.vcf.gz"

STRUCT_IN="/scratch/bjl31194/yaupon/wgs/plates1234/vcf/structure/Ilex384forStructure.recode.strct_in"

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/gwas"

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load modules
ml PLINK/2.0.0-a.6.20-gfbf-2024a
#ml ADMIXTURE/1.3.0
#ml Structure/2.3.4-GCC-12.3.0
#ml Structure_threader/1.3.10-foss-2023a

# move to the proper directory
cd $DATADIR

## Make .bed file for inputting to GEMMA and filter data
## identify prune sites and filter
## KEY: --indep-pairwise x y z
# a) consider a window of x SNPs
# b) calculate LD between each pair of SNPs in the window
# b) remove one of a pair of SNPs if the LD is greater than z
# c) shift the window y SNPs forward and repeat the procedure

plink --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --nonfounders --set-missing-var-ids @:# \
--maf 0.05 --geno 0.1 --mind 0.5 --snps-only \
--make-bed --out Ilex_plates1-5_texas

## attach phenotype data
plink --bfile Ilex_plates1-5_texas --allow-no-sex --pheno texas_sex_phenotypes.txt \
--make-bed --out gemma_input_texas

## Estimating LD with plink

# plink --vcf $VCF --double-id --allow-extra-chr \
# --set-missing-var-ids @:# \
# --maf 0.01 --geno 0.1 --mind 0.5 --chr h1tg000051l \
# -r2 gz --ld-window 100 --ld-window-kb 3000 \
# --ld-window-r2 0 \
# --out Ivom384chr

## identify prune sites
# plink --vcf $VCF --double-id --allow-extra-chr \
# --set-missing-var-ids @:# \
# --indep-pairwise 50 10 0.1 --out Ilex384

# linkage prune and create pca files
# plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
# --extract Ivom384.prune.in \
# --make-bed --pca var-wts --out Ivom384

# # generate structure input file
# plink --bfile Ilex384 --allow-extra-chr --recode structure --out Ilex384forStructure
# # generate  "0,1,2" coded genotype matrix
# plink --bfile Ilex384 --allow-extra-chr --recode A --out Ilex384forRDA

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

#structure_threader run -Klist 2 3 4 5 6 -R 3 -i $STRUCT_IN -o $DATADIR -t 16 --params mainparams_Ivom384 --ind indfile_Ivom384names -st /apps/eb/Structure/2.3.4-GCC-12.3.0/bin/structure
#structure_threader plot -i . -f structure -K 2 3 4 5 6 --ind indfile_Ivom384
## other misc scripts ##

# run pong locally for ADMIXTURE visualization
# use Q matrix files from ADMIXTURE output
#pong -m pong_filemap.txt -i ind2pop.txt
#pong -m pong_filemap_Ilex384.txt -i ind2pop_Ilex384.txt
#structure_threader run -Klist 2 3 4 5 6 -R 3 -i Ivom384forStructure.recode.strct_in -o . -t 16 --params mainparams_Ivom384 --ind indfile.csv -st /apps/eb/Structure/2.3.4-GCC-12.3.0/bin/structure

# structure mainparams:  onerowperind TRUE; label TRUE; popdata, popflag FALSE; locdata FALSE; phenotype FALSE; extracols 1; markernames, mapdistances TRUE