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

Ivom384 <- read.structure("Ivom384forStructureRecode.STR", 
                          n.ind=366,
                          n.loc=36567, 
                          onerowperind=TRUE,
                          col.lab=1,
                          col.pop=0,
                          row.marknames=1,
                          NA.char = "0",) #366 gts, 36567 markers 

summary(Ivom384)

ind2pop <- read.table("ind2pop.txt")
pop_clust <- read.table("pop_clust.txt")
sites <- read.table("sites_Ivom384.txt")
#attach pop info
Ivom384@pop <- as.factor(sites$V1)
Ivom384pop <- genind2genpop(Ivom384)

#calc population allele freqs
pop_freqs <- makefreq(Ivom384pop)
#calc Fst
fst_data <- genind2hierfstat(Ivom384) 
#Nei's (1987) Fst
basic.stats(Ivom384)
wc(Ivom384)
#pairwise Fst
pairwiseFst <- genet.dist(Ivom384, method = "WC84")

# using Kmeans and DAPC in adegenet 
set.seed(29475018) # Setting a seed for a consistent result
clust <- find.clusters(Ivom384, max.n.clust = 10) 
names(grp)
dapc1 <- dapc(Ivom384, clust$grp) 
scatter(dapc1) # plot of the group
biplot(dapc1)

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
        cstar=0, cpoint=3, grid=FALSE)

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




