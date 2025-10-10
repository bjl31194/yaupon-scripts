#get heterozygosity per locus from VCF

bcftools view -Ob Ivom384_no_maf.vcf.gz | bcftools query -f '%CHROM\t%POS[\t%GT]\n' | \
awk '{
  het=0; total=0;
  for (i=3;i<=NF;i++) {
    if ($i=="0/1" || $i=="1/0") het++;
    if ($i!="./.") total++;
  }
  if (total > 0) print $1"\t"$2"\t"het/total;
}' > Ivom384.het