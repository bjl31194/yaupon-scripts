#!/bin/bash
#SBATCH --job-name=bcftools_getcandsnps
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew"

VCF_IN="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ilex1-5_names.vcf.gz"
VCF_OUT="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ilex1-5_vcftools_filter.vcf.gz"

# set filters
MAF=0.01 # mean was 0.03, median 0.005
MISS=0.9 # mean was 0.96 (only missing 4%)
QUAL=30 # 0.01% error rate or lower
MIN_DEPTH=7 # 1st quantile 7.23, mean was 8.3
MAX_DEPTH=60 # went with mean depth x 2 + a little extra judging from where most variants fell on histogram

# load modules
#ml VCFtools/0.1.16-GCC-13.3.0
ml BCFtools/1.21-GCC-13.3.0

# move to the vcf directory
cd $DATADIR

#######################
## using bcftools ##
#######################
# biallelic SNPs with less than 20% missing data, <0.01% error rate, 250X total read depth across samples, and 7-60X depth per sample
# bcftools view -Oz --threads 8 -m2 -M2 -v snps -i 'F_MISSING<0.2 & QUAL > 30 & INFO/DP > 250 & FMT/DP > 7 & FMT/DP < 60' Ilex1-5_names.vcf.gz -o Ilex1-5_names_filter.vcf.gz

# bcftools stats Ilex_plates1-5_merged.vcf.gz > Ilex_merged.stats
# bcftools stats Ilex1-5_pruned.vcf.gz > Ilex_pruned.stats
# bcftools stats Ivom1-5_filter.vcf.gz > Ivom1-5_old.stats

# subset big vcf for mkado
# bcftools view --threads 8 -Oz -S Ivom1-5_newnames.txt Ilex1-5_names_filter.vcf.gz > Ivom1-5_names_newfilter.vcf.gz
# bcftools view --threads 8 -Oz -S outgroup_Ipa.txt Ilex1-5_names_filter.vcf.gz > MC-IP-2.vcf.gz

# bcftools index --threads 8 -t Ivom1-5_names_newfilter.vcf.gz
# bcftools index --threads 8 -t MC-IP-2.vcf.gz

##########################################
## perform filtering with vcftools ##
##########################################
# bcftools reheader --threads 8 --samples ./Ilex1-5_newnames.txt Ilex_plates1-5_merged.vcf.gz -o Ilex1-5_names.vcf.gz

# vcftools --gzvcf $VCF_IN \
# --remove-indels --maf $MAF --minQ $QUAL \
# --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
# --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --max-missing $MISS --recode --stdout | bgzip -c > \
# $VCF_OUT

# optional filtering by sample id to remove decidua individuals (requires txt file with list of sample names to keep)
#bcftools view -Oz -S only_yaupon.txt Ivom_plate1_filter.vcf.gz > Ivom_plate1_sppfilter.vcf.gz

#bcftools stats -S ./namechange_Ivom384.txt $VCF_IN > Ivom384_filtered_stats.vchk
#plot-vcfstats -p . Ivom384_filtered_stats.vchk

# rename samples in vcf header
# bcftools reheader --samples ./Ilex1-5_newnames.txt -o Ilex1-5_names.vcf.gz

## for filtering sexed vomitoria individuals
# bcftools view -Oz -S only_sexed.txt Ilex_plates1-5_names_filter.vcf.gz > Ilex_plates1-5_names_filter_sexed.vcf.gz

## filter out sets of individuals and remove monomorphic sites
# vcftools --gzvcf Ivom1-5_filter.vcf.gz \
# --keep all_florida.txt \
# --recode --recode-INFO-all --stdout | \
# bcftools view -i "MAC>=1" \
# -o Ivom1-5_florida.vcf.gz

# vcftools --gzvcf Ilex_plates1-5_names_filter.vcf.gz \
# --keep redrep.txt \
# --max-missing 1 \
# --recode --recode-INFO-all --stdout | \
# bcftools view -i "MAC>=1" \
# -o Ilex_redrep_nomiss.vcf.gz

# vcftools --gzvcf Ivom1-5_florida.vcf.gz \
# --keep all_florida.txt \
# --window-pi 50000 \
# --out florida_pi_50kb

# ## Tajima's D for genetic subpops
# vcftools --gzvcf Ivom1-5_atlantic.vcf.gz \
# --TajimaD 100000 \
# --out Ivom1-5_atl_TajD_100kb

# vcftools --gzvcf Ivom1-5_gulf.vcf.gz \
# --TajimaD 100000 \
# --out Ivom1-5_gulf_TajD_100kb

# vcftools --gzvcf Ivom1-5_florida.vcf.gz \
# --TajimaD 100000 \
# --out Ivom1-5_fl_TajD_100kb

## query VCF file for specific variants by position
# vcftools --gzvcf Ivom1-5_filter.vcf.gz --chr Chr05 --from-bp 17957200 --to-bp 17957300 --recode --recode-INFO-all --out QTL_Chr05

## query for EHH cand region SNPs
# bcftools view Ivom1-5_filter.vcf.gz -r Chr02:41050000-41540000,Chr03:1170000-1660000,Chr04:6770000-7260000,\
# Chr06:27470000-27960000,Chr12:16970000-17460000,Chr13:6750000-7240000 \
# -Ov -o EHH_cand_region_snps.vcf

#bcftools view Ivom_wild_filter.vcf -Oz -o Ivom_wild_filter.vcf.gz
bcftools index -t --threads 4 Ivom_wild_filter.vcf.gz
bcftools view Ivom_wild_filter.vcf.gz -R cand_regions_EHH_dune.txt \
-Ov -o EHH_dune_cand_region_snps.vcf