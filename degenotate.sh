#!/bin/bash
#SBATCH --job-name=mkado_imputed_Idis
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40gb
#SBATCH --time=3-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error

# set parameters
OUTDIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/mkado"
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew"

REF="/scratch/bjl31194/yaupon/references/Ilex_vomitoria_var_GA_F_4_HAP1_V1_release/Ilex_vomitoria_var_GA_F_4/sequences/Ilex_vomitoria_var_GA_F_4.HAP1.mainGenome.fasta"
GFF="/scratch/bjl31194/yaupon/wgs/plates1-5/ann/Ilex_Hap1.filter.gff3"
VCF="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ivom1-5_names_nofilter.vcf.gz"
OUT="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/MC-IP-02.vcf.gz"
MAF=0.1
GENES="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/gene_IDs.txt"


cd $DATADIR

# index vcf
ml BCFtools/1.21-GCC-13.3.0
ml tabix/0.2.6-GCCcore-13.3.0
# bcftools index -t $VCF --threads 8 

# call variants in outgroup
bcftools mpileup -a AD,DP,SP,INFO/AD -Ou -f $REF /scratch/bjl31194/yaupon/wgs/plates1-5/align_hap1/25055FL-03-01-91_S91_L007.Ivo.sorted.bam | \
    bcftools call --threads 8 -mv -Oz -o ${DATADIR}/MC-IP-02.vcf.gz

# subset Ivom individuals from big unfiltered vcf
# bcftools view --threads 8 -Oz -S Ivom1-5_newnames.txt Ilex1-5_names.vcf.gz > Ivom1-5_names_nofilter.vcf.gz

# index the vcfs
tabix -p vcf MC-IP-02.vcf.gz
# tabix -p vcf /scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/Ivom1-5_names_nofilter.vcf.gz

conda init

# load conda environment
# conda activate degenotate
conda activate mkado

# run degenotate
# /home/bjl31194/.conda/envs/degenotate/bin/degenotate.py -a ${GFF} -g ${REF} -v ${VCF} -maf ${MAF} -e exclude.txt -u outgroup_Ipa.txt -o ${OUTDIR} -sfs --overwrite

# run mkado (asymptotic MK)
# /home/bjl31194/.conda/envs/mkado/bin/mkado vcf --vcf ${VCF} --outgroup-vcf ${OUT} --ref ${REF} --gff ${GFF} --per-gene --workers 8 -f tsv --verbose -a 

# run mkado (regular MK)
/home/bjl31194/.conda/envs/mkado/bin/mkado vcf --vcf ${VCF} --outgroup-vcf ${OUT} --ref ${REF} --gff ${GFF} --per-gene --workers 16 --verbose -f tsv > ${DATADIR}/mkado/mkado_results_Ipa_normal.tsv

## get fasta from candidate regions in GFF format
# ml BEDTools/2.31.1-GCC-13.3.0
# bedtools getfasta -fi /scratch/bjl31194/yaupon/references/Ilex_vomitoria_var_GA_F_4_HAP1_V1_release/Ilex_vomitoria_var_GA_F_4/sequences/Ilex_vomitoria_var_GA_F_4.HAP1.mainGenome.fasta \
# -bed coast_gene_hits.gff3 -fo coast_gene_seqs.fasta
# # split fasta
# awk -v RS='>' -v ORS='>' '(NR%100) == 1 { close(out); out="candidate_gene_seqs"(++n_seq)".fasta" } { print > out }' gulfxatl_candidate_seqs_wg.fasta
