cd /scratch/bjl31194/yaupon/trimmed_reads

for file in /scratch/bjl31194/yaupon/trimmed_reads/*.1.fq_trimmed.fq
do
  d=$(dirname "$file") # get dir name
  echo "${d}"
  f=${file##*/} # remove heading path
  echo "${f}"
  f1=${f%.1.fq_trimmed.fq} # remove trailing
  echo "${f1}"
  sed 's|pqX|'$f1'|g' ~/yaupon/bowtie/bowtie.sh > ~/yaupon/bowtie/subbowtie_$f1.sh  #replace preset keyword ychX
  echo "sbatch ~/yaupon/bowtie/subbowtie_$f1.sh" >> ~/yaupon/bowtie/qbowtie.sh
done
