library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
install.packages("qqman")
library(qqman)

gemmaResults <- read_table("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1234/GWAS_results_sex_noLP.lmm.assoc.txt")       
gemmaResults_texas <- read_table("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1234/GWAS_results_sex_texas.lmm.assoc.txt")
gemmaResults_gulf <- read_table("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1234/GWAS_results_sex_gulf.lmm.assoc.txt")
gemmaResults_atlantic <- read_table("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1234/GWAS_results_sex_atlantic.lmm.assoc.txt")

bonferroni <- -log10(0.05 / nrow(gemmaResults))
bonf_texas <- -log10(0.05 / nrow(gemmaResults_texas))
bonf_gulf <- -log10(0.05 / nrow(gemmaResults_gulf))
bonf_atlantic <- -log10(0.05 / nrow(gemmaResults_atlantic))

manhattan(gemmaResults, chr="chr", bp="ps", p="p_lrt", snp="rs", genomewideline=bonferroni)
manhattan(gemmaResults_texas, chr="chr", bp="ps", p="p_lrt", snp="rs", genomewideline=bonf_texas)
manhattan(gemmaResults_gulf, chr="chr", bp="ps", p="p_lrt", snp="rs", genomewideline=bonf_gulf)
manhattan(gemmaResults_atlantic, chr="chr", bp="ps", p="p_lrt", snp="rs", genomewideline=bonf_atlantic)

top <- gemmaResults %>%
  mutate(negLogP = -log10(p_lrt)) %>%
  select(chr, rs, p_lrt,negLogP) %>%
  filter(negLogP > bonferroni)

top_texas <- gemmaResults_texas %>%
  mutate(negLogP = -log10(p_lrt)) %>%
  select(chr, rs, p_lrt,negLogP) %>%
  filter(negLogP > bonf_texas)

top_gulf <- gemmaResults_gulf %>%
  mutate(negLogP = -log10(p_lrt)) %>%
  select(chr, rs, p_lrt,negLogP) %>%
  filter(negLogP > bonf_gulf)

top_atl <- gemmaResults_atlantic %>%
  mutate(negLogP = -log10(p_lrt)) %>%
  select(chr, rs, p_lrt,negLogP) %>%
  filter(negLogP > bonf_atlantic)

write.csv(top, "yaupon_GWAS_sex_all.csv", row.names = FALSE)
write.csv(top_texas, "yaupon_GWAS_sex_texas.csv", row.names = FALSE)

## plotting snps
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)
sex <- c("females","males")
hom <- c(39,1)
het <- c(1,39)

#atlantic
atlantic <- read.csv("atlantic.csv")

atl <- ggplot(atlantic, aes(call, ..count..)) + 
  geom_bar(aes(fill = sex), position = "stack") +
  coord_cartesian(ylim = c(0, 80)) +
  ggtitle("Atlantic Populations")
atl

#atlantic + FL
atlanticFL <- read.csv("atlanticFL.csv")

atlFL <- ggplot(atlanticFL, aes(call, ..count..)) + 
  geom_bar(aes(fill = sex), position = "stack") +
  coord_cartesian(ylim = c(0, 80)) +
  ggtitle("Atlantic and Florida Pops")
atlFL

#gulf
gulf <- read.csv("gulf.csv")

gf <- ggplot(gulf, aes(call, ..count..)) + 
  geom_bar(aes(fill = sex), position = "stack") +
  coord_cartesian(ylim = c(0, 80)) +
  ggtitle("Gulf Populations")
gf

#all
all <- read.csv("all.csv")

a <- ggplot(all, aes(call, ..count..)) + 
  geom_bar(aes(fill = sex), position = "stack") +
  coord_cartesian(ylim = c(0, 130)) +
  ggtitle("All Populations")
a

#florida
florida <- read.csv("florida.csv")

fl <- ggplot(florida, aes(call, ..count..)) + 
  geom_bar(aes(fill = sex), position = "stack") +
  coord_cartesian(ylim = c(0, 80)) +
  ggtitle("Florida Populations")
fl

