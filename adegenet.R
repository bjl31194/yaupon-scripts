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
library(tidyr)
setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1-5")

Ivom384 <- read.structure("Ivom384forStructureRecode.STR", 
                          n.ind=366,
                          n.loc=36567, 
                          onerowperind=TRUE,
                          col.lab=1,
                          col.pop=0,
                          row.marknames=1,
                          NA.char = "0",) #366 gts, 36567 markers 

IvnoMAF <- read.vcfR("Ivom384_no_maf.vcf.gz", nrows = 92788)
Ivom_vcf <- read.vcfR("Ivom_only_384_filtered_names.vcf", nrows = 92976)
Ivom384_noMAF <- vcfR2genind(IvnoMAF)
ld_pruned_matrix <- read_table("Ivom384.prune.in", col.names = FALSE)
summary(Ivom384)

ind2pop <- read.table("ind2pop.txt")
pop_clust <- read.table("pop_clust.txt")
sites <- read.table("sites_Ivom384.txt")
states <- read.table("states.txt")
#attach pop info
Ivom384@pop <- as.factor(sites$V1)
Ivom384@pop <- as.factor(states$V1)
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

## calculate latlong correlation with first two PCs ##
pca_latlong <- coords %>% 
  mutate(PC1 = pca$PC1) %>%
  mutate(PC2 = pca$PC2)

cor(pca_latlong)

## Pairwise Fst heatmap ##
require(ape)
pwFst <- read.csv("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1234/pairwiseFst_Ivom384.csv", row.names = 1)

pwFst2 <- pwFst %>%
  rownames_to_column() %>%
  gather(colname, Fst, -rowname)

pwFst2$colname <- gsub('\\.', '-', pwFst2$colname)

colnames(pwFst2) <- c("pop1","pop2","Fst")

pwFst_heat <- ggplot(pwFst2, aes(x = pop1, y = pop2, fill = Fst)) +
  geom_tile()
pwFst_heat

pwFst <- as.matrix(pwFst)

Ivom384.tree <- nj(pwFst)
plot(Ivom384.tree, type="unr", tip.col=funky(nPop(Ivom384)), font=2, )
annot <- round(Ivom384.tree$edge.length,2)
edgelabels(annot[annot>0], which(annot>0), frame="n")
add.scale.bar()

table.paint(pwFst, col.labels=colnames(pwFst))

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

## nucleotide diversity

gulf_pi <- read_table("gulf_pi_50kb.windowed.pi")
atl_pi <- read_table("atlantic_pi_50kb.windowed.pi")
fl_pi <- read_table("florida_pi_50kb.windowed.pi")

gulf_pi <- gulf_pi %>%
  mutate(CHROM = gsub("[a-zA-Z ]", "", gulf_pi$CHROM)) 
atl_pi <- atl_pi %>%
  mutate(CHROM = gsub("[a-zA-Z ]", "", atl_pi$CHROM))
fl_pi <- fl_pi %>%
  mutate(CHROM = gsub("[a-zA-Z ]", "", fl_pi$CHROM))

gulf_pi$CHROM <- as.numeric(gulf_pi$CHROM)
atl_pi$CHROM <- as.numeric(atl_pi$CHROM)
fl_pi$CHROM <- as.numeric(fl_pi$CHROM)

gulf_avg_pi <- mean(gulf_pi$PI)
atl_avg_pi <- mean(atl_pi$PI)
fl_avg_pi <- mean(fl_pi$PI)

manhattan(gulf_pi, chr="CHROM", snp="N_VARIANTS", bp="BIN_START", p="PI", logp=FALSE, ylim=c(0,0.0015))
manhattan(atl_pi, chr="CHROM", snp="N_VARIANTS", bp="BIN_START", p="PI", logp=FALSE, ylim=c(0,0.0015))
manhattan(fl_pi, chr="CHROM", snp="N_VARIANTS", bp="BIN_START", p="PI", logp=FALSE, ylim=c(0,0.0015))

## Tajima's D ##
require(qqman)
require(tidyr)

Ivom_TajD <- read_table("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1-5/Ivom1-5_TajD_50kb.Tajima.D")

Ivom_TajD <- Ivom_TajD %>%
  mutate(CHROM = gsub("[a-zA-Z ]", "", Ivom_TajD$CHROM))

Ivom_TajD$CHROM <- as.numeric(Ivom_TajD$CHROM)

Ivom_TajD$TajimaD <- na_if(Ivom_TajD$TajimaD, "nan")

Ivom_TajD <- Ivom_TajD %>%
  drop_na(TajimaD)

Ivom_TajD$TajimaD <- as.numeric(Ivom_TajD$TajimaD)

Ivom_TajD <- Ivom_TajD %>%
  filter(N_SNPS > 20)

manhattan(Ivom_TajD, chr="CHROM", snp="N_SNPS", bp="BIN_START", p="TajimaD", logp=FALSE, ylim=c(-5,8))

Ivom_avg_D <- mean(Ivom_TajD$TajimaD)


ggplot(data=Ivom_TajD)  +
         geom_line(aes(x=N_SNPS,y=TajimaD))

# by pop
gulf_D <- read_table("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1-5/Ivom1-5_gulf_TajD_100kb.Tajima.D")
atl_D <- read_table("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1-5/Ivom1-5_atl_TajD_100kb.Tajima.D")
fl_D <- read_table("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1-5/Ivom1-5_fl_TajD_100kb.Tajima.D")

gulf_D <- gulf_D %>%
  mutate(CHROM = gsub("[a-zA-Z ]", "", gulf_D$CHROM)) 
atl_D <- atl_D %>%
  mutate(CHROM = gsub("[a-zA-Z ]", "", atl_D$CHROM))
fl_D <- fl_D %>%
  mutate(CHROM = gsub("[a-zA-Z ]", "", fl_D$CHROM))

gulf_D$CHROM <- as.numeric(gulf_D$CHROM)
atl_D$CHROM <- as.numeric(atl_D$CHROM)
fl_D$CHROM <- as.numeric(fl_D$CHROM)

gulf_D$TajimaD <- na_if(gulf_D$TajimaD, "nan")
gulf_D <- gulf_D %>%
  drop_na(TajimaD)
gulf_D$TajimaD <- as.numeric(gulf_D$TajimaD)

atl_D$TajimaD <- na_if(atl_D$TajimaD, "nan")
atl_D <- atl_D %>%
  drop_na(TajimaD)
atl_D$TajimaD <- as.numeric(atl_D$TajimaD)

fl_D$TajimaD <- na_if(fl_D$TajimaD, "nan")
fl_D <- fl_D %>%
  drop_na(TajimaD)
fl_D$TajimaD <- as.numeric(fl_D$TajimaD)

DATASET <- fl_D

manhattan(DATASET, chr="CHROM", snp="N_SNPS", bp="BIN_START", p="TajimaD", logp=FALSE, ylim=c(-5,8))

gulf_avg_D <- mean(gulf_D$TajimaD)
atl_avg_D <- mean(atl_D$TajimaD)
fl_avg_D <- mean(fl_D$TajimaD)

ggplot(data=Ivom_TajD)  +
  geom_line(aes(x=N_SNPS,y=TajimaD))

# NJ tree using ape
redrepvcf <- read.vcfR("Ilex_redrep_filter.vcf.gz")
redrepbin <- vcfR2DNAbin(redrepvcf, extract.haps = FALSE, unphased_as_NA = FALSE, verbose = TRUE)
dist <- dist.dna(redrepbin, model = "TN93", gamma = TRUE)
nj <- njs(dist)
plot.phylo(nj, type = "phylogram", use.edge.length = FALSE, cex=0.5)


ape::image.DNAbin(Ilex_redrep_align[,ape::seg.sites(Ilex_redrep_align)],cex.lab=0.5)

heatmap <- as.data.frame(as.matrix(dist))
table.paint(heatmap, cleg=0, clabel.row=.5, clabel.col=.5)


##########
## SNMF ##
##########

install.packages("BiocManager")
BiocManager::install("LEA", force = TRUE)
library(LEA)

### BEGIN Kelly Petersen Code ###

genind2structure <- function(obj, file="", pops=TRUE){
if(!"genind" %in% class(obj)){
  warning("Function was designed for genind objects.")
}

# get the max ploidy of the dataset
pl <- max(obj@ploidy)
# get the number of individuals
S <- adegenet::nInd(obj)
# column of individual names to write; set up data.frame
tab <- data.frame(ind=rep(adegenet::indNames(obj), each=pl))
# column of pop ids to write
if(pops){
  popnums <- 1:adegenet::nPop(obj)
  names(popnums) <- as.character(unique(adegenet::pop(obj)))
  popcol <- rep(popnums[as.character(adegenet::pop(obj))], each=pl)
  tab <- cbind(tab, data.frame(pop=popcol))
}
loci <- adegenet::locNames(obj) 
# add columns for genotypes
tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=adegenet::nLoc(obj),
                         dimnames=list(NULL,loci)))

# begin going through loci
for(L in loci){
  thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), 
                            dimnames(obj@tab)[[2]]), 
                      drop = FALSE] # genotypes by locus
  al <- 1:dim(thesegen)[2] # numbered alleles
  for(s in 1:S){
    if(all(!is.na(thesegen[s,]))){
      tabrows <- (1:dim(tab)[1])[tab[[1]] == adegenet::indNames(obj)[s]] # index of rows in output to write to
      tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
      tab[tabrows,L] <- rep(al, times = thesegen[s,])
    }
  }
}

# export table
write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
}

genind2structure(Ivom384, file="structure_Ivom384.txt", pops=TRUE)



### END Kelly Petersen Code ###

geno <- read.geno("Ivom384.geno")

project = snmf("Ivom384.geno",
               K = 1:10, 
               entropy = TRUE, 
               repetitions = 10,
               project = "new")

# plot cross-entropy criterion of all runs of the project
plot(project, cex = 1.2, col = "lightblue", pch = 19)

# show the project
show(project)

# summary of the project
summary(project)

# get the cross-entropy of all runs for K = 4
ce = cross.entropy(project, K = 4)

# select the run with the lowest cross-entropy for K = 4
best = which.min(ce)

# display the Q-matrix
Q.matrix <- as.qmatrix(Q(project, K = 4, run = best))
my.colors <- c("tomato", "lightblue", "olivedrab", "gold")

barplot(Q.matrix, 
        border = NA, 
        space = 0, 
        col = my.colors, 
        xlab = "Individuals",
        ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp

axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .4)


# get the ancestral genotype frequency matrix, G, for the 2nd run for K = 4. 
G.matrix = G(project, K = 4, run = 2)

###
data("tutorial")
write.geno(tutorial.R, "genotypes.geno")


tutorial = snmf("genotypes.geno",
               K = 1:10, 
               entropy = TRUE, 
               repetitions = 10,
               project = "new")

# plot cross-entropy criterion of all runs of the project
plot(tutorial, cex = 1.2, col = "lightblue", pch = 19)

## plot He by Ho

temp <- summary(Ivom384)

plot(temp$Hexp, temp$Hobs, pch=5, cex=1, xlim=c(0,1), ylim=c(0,1))
abline(0,1,lty=2)
