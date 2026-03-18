#!/bin/bash
#SBATCH --job-name=plink_Ilex
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=3-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew"

VCF="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ilex1-5_names_filter.vcf.gz"

PREFIX="Ilex1-5"

SUBSET="atlantic"

STRUCT_IN="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/structure/Ivom1-5_forStructure.recode.strct_in"

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/"

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load modules
ml PLINK/2.0.0-a.6.20-gfbf-2024a
#ml ADMIXTURE/1.3.0
#ml Structure/2.3.4-GCC-12.3.0
#ml Structure_threader/1.3.10-foss-2023a
#ml Beagle/5.4.22Jul22.46e-Java-11

# move to the proper directory
cd $OUTDIR

## Make .bed file for inputting to GEMMA and filter data

# plink --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --nonfounders --set-missing-var-ids @:# \
# --maf 0.0013 --geno 0.1 --mind 0.5 --snps-only \
# --make-bed --out Ilex_plates1-5_${SUBSET}

# ## attach phenotype data
# plink --bfile Ilex_plates1-5_${SUBSET} --allow-no-sex --pheno ${SUBSET}_sex_phenotypes.txt \
# --make-bed --out gemma_input_${SUBSET}

## statistical phasing with BEAGLE on Sapelo2 cluster:

#java -jar ${EBROOTBEAGLE}/beagle.jar gt=Ivom_Ipa_outgroup.vcf.gz nthreads=8 out=Ivom_Ipa_outgroup_phased
#java -jar ${EBROOTBEAGLE}/beagle.jar gt=Ivom1-5_inland.vcf.gz nthreads=8 out=Ivom1-5_inland_phased 

## get fasta from candidate regions in GFF format
# ml BEDTools/2.31.1-GCC-13.3.0
# bedtools getfasta -fi /scratch/bjl31194/yaupon/references/Ilex_vomitoria_var_GA_F_4_HAP1_V1_release/Ilex_vomitoria_var_GA_F_4/sequences/Ilex_vomitoria_var_GA_F_4.HAP1.mainGenome.fasta \
# -bed coast_gene_hits.gff3 -fo coast_gene_seqs.fasta
# # split fasta
# awk -v RS='>' -v ORS='>' '(NR%100) == 1 { close(out); out="candidate_gene_seqs"(++n_seq)".fasta" } { print > out }' gulfxatl_candidate_seqs_wg.fasta

## Estimating LD with plink

# plink --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --nonfounders \
# --set-missing-var-ids @:# \
# --maf 0.01 --geno 0.2 --mind 0.5 --chr Chr02 \
# --thin 0.5 -r2 gz --ld-window 100 --ld-window-kb 1000 \
# --ld-window-r2 0 \
# --make-bed --out Ivom_chr2

## identify prune sites, LD prune, filter variants, and create bed, pca, and structure files
## KEY: --indep-pairwise x y z
# a) consider SNPs in a window of x kb
# b) calculate LD between each pair of SNPs in the window
# b) remove one of a pair of SNPs if the LD is greater than z
# c) shift the window y SNPs forward and repeat the procedure

plink --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
--maf 0.0013 --geno 0.1 --snps-only \
--make-bed --out $PREFIX

plink -bfile $PREFIX --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
--indep-pairwise 100 5 0.7 \
--out $PREFIX

plink --vcf $VCF --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
--extract ${PREFIX}.prune.in \
--make-bed --pca var-wts --out ${PREFIX}_pruned

# plink --bfile ${PREFIX}_pruned --allow-extra-chr --allow-no-sex --recode structure --out ${PREFIX}_forStructure

## generate  "0,1,2" coded genotype matrix
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

# structure_threader run -Klist 2 3 4 5 6 7 8 -R 3 -i $STRUCT_IN -o $OUTDIR -t 16 --params mainparams_Ivom1-5 --ind indfile_Ivom1-5.csv -st /apps/eb/Structure/2.3.4-GCC-12.3.0/bin/structure
# structure_threader plot -i . -f structure -K 2 3 4 5 6 7 8 --ind indfile_Ivom1-5.csv
 
## other misc scripts ##

# run pong locally for ADMIXTURE visualization
# use Q matrix files from ADMIXTURE output
#pong -m pong_filemap.txt -i ind2pop.txt
# pong -m pong_filemap_Ivom1-5_K4.txt -i ind2pop_Ivom1-5.txt -n pop_order_Ivom1-5.txt
#structure_threader run -Klist 2 3 4 5 6 -R 3 -i Ivom384forStructure.recode.strct_in -o . -t 16 --params mainparams_Ivom384 --ind indfile.csv -st /apps/eb/Structure/2.3.4-GCC-12.3.0/bin/structure

# structure mainparams:  onerowperind TRUE; label TRUE; popdata, popflag FALSE; locdata FALSE; phenotype FALSE; extracols 1; markernames, mapdistances TRUE