library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
install.packages("qqman")
library(qqman)

gemmaResults <- read_table("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1234/GWAS_results_sex.lmm.assoc.txt")

bonferroni <- -log10(0.05 / nrow(gemmaResults))

manhattan(gemmaResults, chr="chr", bp="ps", p="p_lrt", snp="rs", genomewideline=bonferroni)

top <- gemmaResults %>%
  mutate(negLogP = -log10(p_lrt)) %>%
  select(chr, rs, p_lrt,negLogP) %>%
  filter(negLogP > bonferroni)
