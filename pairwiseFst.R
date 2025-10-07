library("adegenet")
library("hierfstat")

setwd("/scratch/bjl31194/yaupon/wgs/plates1234/vcf")

Ivom384 <- read.structure("Ivom384forStructureRecode.STR", 
                          n.ind=366,
                          n.loc=36567, 
                          onerowperind=TRUE,
                          col.lab=1,
                          col.pop=0,
                          row.marknames=1,
                          NA.char = "0",) #366 gts, 36567 markers 
sites <- read.table("sites_Ivom384.txt")

Ivom384@pop <- as.factor(sites$V1)

#pairwise Fst
pairwiseFst <- genet.dist(Ivom384, method = "WC84")

write.csv(pairwiseFst, "./pairwiseFst_Ivom384.csv")
