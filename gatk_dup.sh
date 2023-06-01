cd /scratch/bjl31194/yaupon/trimmed_reads

for file in /scratch/bjl31194/yaupon/trimmed_reads/*.Gr.sorted.bam
do
  d=$(dirname "$file") # get dir name
  echo "${d}"
  f=${file##*/} # remove heading path
  echo "${f}"
  f1=${f%.Gr.sorted.bam} # remove trailing
  echo "${f1}"
  sed 's|Ivo|'$f1'|g' ~/yaupon/gatk/gatk_hapcaller.sh > ~/yaupon/gatk/sub_gatk_hapcaller_$f1.sh  #replace preset keyword
  echo "sbatch ~/yaupon/gatk/sub_gatk_hapcaller_$f1.sh" >> ~/yaupon/gatk/q_gatk_hapcaller.sh
done
