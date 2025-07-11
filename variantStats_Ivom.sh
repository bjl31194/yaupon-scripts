#!/bin/bash
#SBATCH --job-name=variantStats_Ilex
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1234/vcf/vcftools"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

DATADIR="/scratch/bjl31194/yaupon/wgs/plates1234/vcf"

# load modules
ml BCFtools/1.18-GCC-12.3.0
ml vcflib/1.0.9-gfbf-2024a-R-4.4.2
ml VCFtools/0.1.16-GCC-13.3.0

cd $DATADIR

# randomly sample VCF (target ~ variants)
bcftools view Ilex_plates1234_merged.vcf.gz | vcfrandomsample -r 0.005 > Ilex_plates1234_subset.vcf

### bcftools plot-vcfstats (does not work) ###

# generate stats
#bcftools stats Ilex_plates1234_subset.vcf > Ilex_plates1234_subset.vchk

# plot stats
#plot-vcfstats -p ./stats Ilex_plates1234_subset.vchk

###

### using vcftools (works) ###

# compress vcf
bgzip Ilex_plates1234_subset.vcf
# index vcf
bcftools index Ilex_plates1234_subset.vcf.gz

SUBSET='/scratch/bjl31194/yaupon/wgs/plates1234/vcf/Ilex_plates1234_subset.vcf.gz'

# calculate allele freqs
vcftools --gzvcf $SUBSET --freq2 --out $OUTDIR --max-alleles 2

# mean depth of coverage
vcftools --gzvcf $SUBSET --depth --out $OUTDIR

# mean depth per site
vcftools --gzvcf $SUBSET --site-mean-depth --out $OUTDIR

# per site qual scores
vcftools --gzvcf $SUBSET --site-quality --out $OUTDIR

# missing data per sample
vcftools --gzvcf $SUBSET --missing-indv --out $OUTDIR

# missing data per site
vcftools --gzvcf $SUBSET --missing-site --out $OUTDIR

# heterozygosity per sample
vcftools --gzvcf $SUBSET --het --out $OUTDIR

