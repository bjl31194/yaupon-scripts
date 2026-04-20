#!/bin/bash
#SBATCH --job-name=subset_beagle
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=3-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# name=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/bjl31194/yaupon/wgs/plate1/read_array.txt)

# command for making read array file:
# ls -1 | sed 's/_L006_R.*//' | uniq > read_array.txt

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew"

VCF="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ivom1-5_names_nofilter.vcf.gz"

PREFIX="Ilex1-5"

SUBSET="atlantic"

STRUCT_IN="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/structure/Ivom1-5_forStructure.recode.strct_in"

OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/"

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load modules
ml Beagle/5.4.22Jul22.46e-Java-11
ml BCFtools/1.21-GCC-13.3.0

# move to the proper directory
cd $OUTDIR

## subset populations from VCF ##
bcftools view -Oz -S inland_atl.txt Ivom_wild_filter.vcf.gz > inland_atl.vcf.gz
bcftools view -Oz -S dune_atl.txt Ivom_wild_filter.vcf.gz > dune_atl.vcf.gz
bcftools view -Oz -S inland_gulf.txt Ivom_wild_filter.vcf.gz > inland_gulf.vcf.gz
bcftools view -Oz -S dune_gulf.txt Ivom_wild_filter.vcf.gz > dune_gulf.vcf.gz


## statistical phasing with BEAGLE ##
java -jar ${EBROOTBEAGLE}/beagle.jar gt=inland_atl.vcf.gz nthreads=8 out=inland_atl_phased
java -jar ${EBROOTBEAGLE}/beagle.jar gt=dune_atl.vcf.gz nthreads=8 out=dune_atl_phased 
java -jar ${EBROOTBEAGLE}/beagle.jar gt=inland_gulf.vcf.gz nthreads=8 out=inland_gulf_phased
java -jar ${EBROOTBEAGLE}/beagle.jar gt=dune_gulf.vcf.gz nthreads=8 out=dune_gulf_phased 
