install.packages("adegenet", dep=TRUE)
library(devtools)
install_github("jgx65/hierfstat")
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library("hierfstat")

setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1234")

Ivom384 <- read.structure("Ivom384forStructureRecode.STR")
summary(Ivom384)

df <- read.table("ind2pop.txt")

Ivom384@pop <- as.factor(df$V1)
Ivom384.hwt <- hw.test(Ivom384, B=0)
Ivom384.hwt

Fst_by_loci <- Fst(as.loci(Ivom384))
