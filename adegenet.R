install.packages("adegenet", dep=TRUE)
library(devtools)
devtools::install_github("jgx65/hierfstat")
devtools::install_github("pievos101/PopGenome")
install.packages("vcfR")
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library("hierfstat")
library(vcfR)

setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1234")

Ivom384 <- read.structure("Ivom384forStructureRecode.STR", 
                          n.ind=366,
                          n.loc=36567, 
                          onerowperind=TRUE,
                          col.lab=1,
                          col.pop=0,
                          row.marknames=1,
                          NA.char = "0",) #366 gts, 36567 markers 

IvnoMAF <- read.vcfR("Ivom384_no_maf.vcf.gz", nrows = 92788)
Ivom384_noMAF <- vcfR2genind(IvnoMAF)
ld_pruned_matrix <- read_table("Ivom384.prune.in", col.names = FALSE)
summary(Ivom384)

ind2pop <- read.table("ind2pop.txt")
pop_clust <- read.table("pop_clust.txt")
sites <- read.table("sites_Ivom384.txt")
states <- read.table("states.txt")
#attach pop info
Ivom384@pop <- as.factor(sites$V1)
Ivom384_noMAF@pop <- as.factor(sites$V1)
strata(Ivom384_noMAF) <- data.frame(Ivom384_noMAF@pop, pop_clust$V1, states$V1)
nameStrata(Ivom384_noMAF) <- ~pop/clust/state
head(strata(Ivom384_noMAF, ~pop/clust/state, combine = FALSE))
Ivom384pop <- genind2genpop(Ivom384)
Ivom384pop_noMAF <- genind2genpop(Ivom384_noMAF)
fivepops <- Ivom384[Ivom384@pop %in% c("AR-1","AR-3","VA-1","NC-1")]

#calc population allele freqs
pop_freqs <- makefreq(Ivom384pop)
pop_freqs_noMAF <- makefreq(Ivom384pop_noMAF)
#calc Fst
fst_data <- genind2hierfstat(Ivom384) 
#Nei's (1987) Fst
basic.stats(Ivom384)
wc(Ivom384)
#pairwise Fst
pairwiseFst <- genet.dist(Ivom384, method = "WC84")
test_Fst <- genet.dist(fivepops, method = "WC84")
testFst_matrix <- as.matrix(test_Fst)
write.csv(testFst_matrix, "./testFst.csv")

## PCA ##
# LD prune
Ivom384_pruned <- Ivom384[, loc = ld_pruned_matrix$V1]
# remove outliers (optional)
Ivom384_no_outliers <- Ivom384_pruned[, loc = no_outliers$loci]

sum(is.na(Ivom384_pruned$tab))
obj <- missingno(Ivom384_pruned, type = "mean", cutoff = 0.20, quiet = FALSE, freq = FALSE)
obj2 <- missingno(Ivom384_no_outliers, type = "mean", cutoff = 0.20, quiet = FALSE, freq = FALSE)
pca1 <- dudi.pca(obj$tab, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 3)
pca2 <- dudi.pca(obj2$tab, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:50], main = "Eigenvalues")

s.class(pca1$li, obj$pop, lab = obj$pops, sub = "PCA 1-2", csub = 2)
add.scatter.eig(pca1$eig[1:20], nf = 3, xax = 1, yax = 2, posi = "top")

s.class(pca1$li, obj$pop, xax = 1, yax = 3, lab = obj$pop, sub = "PCA 1-3", csub = 2)

# using Kmeans and DAPC in adegenet 
set.seed(29475018) # Setting a seed for a consistent result
clust <- find.clusters(Ivom384_pruned, max.n.clust = 10) 
names(grp)
dapc1 <- dapc(Ivom384_pruned, clust$grp) 
scatter(dapc1) # plot of the group

contrib <- loadingplot(dapc1$var.contr, axis=2, thres=.07, lab.jitter=1)


Ivom384.hwt <- hw.test(Ivom384, B=0)
Ivom384.hwt

Fst_by_loci <- Fst(as.loci(Ivom384))


X <- scaleGen(Ivom384, NA.method="mean")

pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

s.class(pca1$li, pop(Ivom384))
title("PCA of Ivom dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li,pop(Ivom384),xax=1,yax=3,sub="PCA 1-3",csub=2)
title("PCA of microbov dataset\naxes 1-3")
add.scatter.eig(pca1$eig[1:20],nf=3,xax=1,yax=3)

col <- funky(15)
s.class(pca1$li, pop(Ivom384),xax=1,yax=3, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE, clabel = 0.5)

s.class(pca1$li, pop(Ivom384),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE, clabel = 0.4)

colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA of microbov dataset\naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

grp <- find.clusters(Ivom384, max.n.clust=20)
dapc1 <- dapc(Ivom384, grp$grp)
scatter(dapc1)
contrib <- loadingplot(dapc1$var.contr, axis=2,
                       thres=.07, lab.jitter=1)
assignplot(dapc1, subset=1:49)

gen <- as(Ivom384, "matrix")
sum(is.na(gen))

## calculate nucleotide diversity (pi)

genind2loci(Ivom384)

## calculate latlong correlation with first two PCs ##
pca_latlong <- coords %>% 
  mutate(PC1 = pca$PC1) %>%
  mutate(PC2 = pca$PC2)

cor(pca_latlong)

## Pairwise Fst heatmap ##
pwFst <- read.csv("pairwiseFst_Ivom384.csv", row.names = 1)

pwFst2 <- pwFst %>%
  rownames_to_column() %>%
  gather(colname, Fst, -rowname)

pwFst2$colname <- gsub('\\.', '-', pwFst2$colname)

colnames(pwFst2) <- c("pop1","pop2","Fst")

pwFst_heat <- ggplot(pwFst2, aes(x = pop1, y = pop2, fill = Fst)) +
  geom_tile()
pwFst_heat

## count private alleles ##
install.packages("poppr")
library(poppr)

pralleles <- private_alleles(Ivom384_noMAF, form = alleles ~ clust, report = "data.frame")
private_alleles(Ivom384_noMAF, form = alleles ~ clust, report = "data.frame")

data(Pinf)
private_alleles(Pinf)
(pal <- private_alleles(Pinf, locus ~ Country, count.alleles = FALSE))
sweep(pal, 2, nAll(Pinf)[colnames(pal)], FUN = "/")
Pinfpriv <- private_alleles(Pinf, report = "data.frame")
ggplot(pralleles) + geom_tile(aes(x = population, y = allele, fill = count))

private_alleles_by_state <- private_alleles(Ivom384_noMAF, form = alleles ~ state, report = "data.frame")

## SNMF ##

install.packages("BiocManager")
BiocManager::install("LEA", force = TRUE)
library(LEA)

vcf2geno("Ivom384_filtered_names.vcf", "Ivom384")

