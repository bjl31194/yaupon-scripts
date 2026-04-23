#!/bin/bash
#SBATCH --job-name=structure_wild
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=7-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

## set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/structure/wild_only"

VCF="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ivom_wild_filter.vcf.gz"

PREFIX="Ivom1-5_wild"

STRUCT_IN="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/structure/Ivom_wild_forStructure.recode.strct_in"

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/structure/wild_only"

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

## load modules
ml PLINK/2.0.0-a.6.20-gfbf-2024a
ml Structure/2.3.4-GCC-12.3.0
ml Structure_threader/1.3.10-foss-2023a

## move to the proper directory
cd $OUTDIR

## PLINK scripts for generating structure input file from VCF

## This line performs basic site filtering for missing data, as well as linkage pruning to remove redundant SNPs that are highly linked. 
## The various options are mostly to help plink deal with non-human data
## How the linkage pruning works: --indep-pairwise x y z
# a) consider a window of x SNPs
# b) calculate LD between each pair of SNPs in the window
# b) remove one of a pair of SNPs if the LD is greater than z
# c) shift the window y SNPs forward and repeat the procedure
plink --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5 --snps-only \
--indep-pairwise 50 10 0.5 \
--out $PREFIX

## This line uses the list of pruned/filtered SNPs from the previous line and makes a BED file from them, as well as eigenvalues/vectors for PCA
plink --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
--extract ${PREFIX}.prune.in \
--make-bed --pca var-wts --out $PREFIX

## This takes the BED file and spits out a (mostly) structure-formatted file 
plink --bfile $PREFIX --allow-extra-chr --allow-no-sex --recode structure --out ${PREFIX}_forStructure


## STRUCTURE - for running on cluster ##
## structure mainparams:  onerowperind TRUE; label TRUE; popdata, popflag FALSE; locdata FALSE; \
## phenotype FALSE; extracols 1; markernames, mapdistances TRUE
structure_threader run -Klist 2 3 4 5 6 7 8 9 10 -R 3 -i $STRUCT_IN -o $OUTDIR -t 16 --params mainparams_Ivom_wild --ind indfile_Ivom_wild.csv -st /apps/eb/Structure/2.3.4-GCC-12.3.0/bin/structure
structure_threader plot -i . -f structure -K 2 3 4 5 6 7 8 9 10 --ind indfile_Ivom_wild.csv

## pong - can install via conda and use locally for STRUCTURE/ADMIXTURE visualization
## use the admixture proportions from STRUCTURE output (*.f)
# pong -m pong_filemap_Ivom1-5_K4.txt -i ind2pop_Ivom1-5.txt -n pop_order_Ivom1-5.txt