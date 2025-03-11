install.packages("vcfR")
library(vcfR)
install.packages("vegan")
library(vegan)
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
library(ggrepel)
setwd("~/yaupon/vcf")

snps <- vcfR::read.vcfR("cohort1_biallelic_QD10_SNPs.vcf", convertNA = TRUE)

snps_num <- vcfR::extract.gt(snps, 
                             element = "GT",
                             IDtoRowNames  = F,
                             as.numeric = T,
                             convertNA = T,
                             return.alleles = F)

snps_num_t <- t(snps_num) 
snps_num_df <- data.frame(snps_num_t) 

snps_num_t <- snps_num_t[-c(7,25,30,35),]

find_NAs <- function(x){
  NAs_TF <- is.na(x)
  i_NA <- which(NAs_TF == TRUE)
  N_NA <- length(i_NA)
  
  cat("Results:",N_NA, "NAs present\n.")
  return(i_NA)
}

# N_rows
# number of rows (individuals)
N_rows <- nrow(snps_num_t)

# N_NA
# vector to hold output (number of NAs)
N_NA   <- rep(x = 0, times = N_rows)

# N_SNPs
# total number of columns (SNPs)
N_SNPs <- ncol(snps_num)

# the for() loop
for(i in 1:N_rows){
  
  # for each row, find the location of
  ## NAs with snps_num_t()
  i_NA <- find_NAs(snps_num[i,]) 
  
  # then determine how many NAs
  ## with length()
  N_NA_i <- length(i_NA)
  
  # then save the output to 
  ## our storage vector
  N_NA[i] <- N_NA_i
}

# 50% of N_SNPs
cutoff50 <- N_SNPs*0.5

hist(N_NA)            
abline(v = cutoff50, 
       col = 2, 
       lwd = 2, 
       lty = 2)

percent_NA <- N_NA/N_SNPs*100

# Call which() on percent_NA
i_NA_50percent <- which(percent_NA < 50) 

snps_num_t02 <- snps_num_t[i_NA_50percent, ]

#removes invariant columns
invar_omit <- function(x){
  cat("Dataframe of dim",dim(x), "processed...\n")
  sds <- apply(x, 2, sd, na.rm = TRUE)
  i_var0 <- which(sds == 0)
  
  
  cat(length(i_var0),"columns removed\n")
  
  if(length(i_var0) > 0){
    x <- x[, -i_var0]
  }
  
  ## add return()  with x in it
  return(x)                      
}

snps_no_invar <- invar_omit(snps_num_t02) 

#no NAs (replace w mean = imputing)
snps_noNAs <- snps_no_invar

N_col <- ncol(snps_no_invar)
for(i in 1:N_col){
  
  # get the current column
  column_i <- snps_noNAs[, i]
  
  # get the mean of the current column
  mean_i <- mean(column_i, na.rm = TRUE)
  
  # get the NAs in the current column
  NAs_i <- which(is.na(column_i))
  
  # record the number of NAs
  N_NAs <- length(NAs_i)
  
  # replace the NAs in the current column
  column_i[NAs_i] <- mean_i
  
  # replace the original column with the
  ## updated columns
  snps_noNAs[, i] <- column_i
  
}

write.csv(snps_noNAs, file = "SNPs_cleaned.csv",
          row.names = F)

### PCA ###

SNPs_cleaned <- read.csv(file = "SNPs_cleaned.csv")

SNPs_scaled <- scale(SNPs_cleaned)

pca_scaled <- prcomp(SNPs_scaled)

screeplot(pca_scaled, 
          ylab  = "Relative importance",
          main = "SNP data screeplot")

summary_out_scaled <- summary(pca_scaled)

PCA_variation <- function(pca_summary, PCs = 2){
  var_explained <- pca_summary$importance[2,1:PCs]*100
  var_explained <- round(var_explained,1)
  return(var_explained)
}

var_out <- PCA_variation(summary_out_scaled,PCs = 10)

N_columns <- ncol(SNPs_scaled)
barplot(var_out,
        main = "Percent variation Scree plot",
        ylab = "Percent variation explained")
abline(h = 1/N_columns*100, col = 2, lwd = 2)

biplot(pca_scaled)

pca_scores <- vegan::scores(pca_scaled)

pop_id <- c("AR", "AR", "VA", "NC", 
            "AR", "AR", "AR", "VA", "NC",
            "AR", "AR", "VA", "VA", "NC",
            "AR", "AR", "VA", "NC", "NC",
            "AR", "AR", "VA", "NC", "NC",
            "AR", "AR", "VA", "NC", "NC",
            "AR", "AR", "VA", "NC", "NC",
            "AR", "AR", "VA", "NC", "NC")

pop_id_allIV <- c("AR", "AR", "VA", "NC", 
                  "AR", "AR", "VA", "NC",
                  "AR", "AR", "VA", "VA", "NC",
                  "AR", "AR", "VA", "NC", "NC",
                  "AR", "AR", "VA", "NC", "NC",
                  "AR", "VA", "NC", "NC",
                  "AR", "VA", "NC", "NC",
                  "AR", "VA", "NC", "NC")

ind_id <- c("AR11", "AR38", "VA21", "NC23", 
            "AR12", "AR31", "AR41", "VA22", "NC24",
            "AR13", "AR32", "VA11", "VA23", "NC25",
            "AR14", "AR33", "VA12", "NC11", "NC26",
            "AR15", "AR34", "VA13", "NC12", "NC27",
            "AR21", "AR35", "VA14", "NC13", "NC28",
            "AR22", "AR36", "VA15", "NC21", "NC31",
            "AR23", "AR37", "VA16", "NC22", "NC32")

ind_id_allIV <- c("AR11", "AR38", "VA21", "NC23", 
            "AR12", "AR31", "VA22", "NC24",
            "AR13", "AR32", "VA11", "VA23", "NC25",
            "AR14", "AR33", "VA12", "NC11", "NC26",
            "AR15", "AR34", "VA13", "NC12", "NC27",
            "AR35", "VA14", "NC13", "NC28",
            "AR36", "VA15", "NC21", "NC31",
            "AR37", "VA16", "NC22", "NC32")

pca_scores2 <- data.frame(pop_id_allIV, pca_scores)            

ggpubr::ggscatter(data = pca_scores2,
                  y = "PC2",
                  x = "PC1",
                  #color = "pop_id_allIV",
                  #shape = "pop_id_allIV",
                  xlab = "PC1 (9.5% variation)",
                  ylab = "PC2 (6.3% variation)",
                  main = "Yaupon SNPs PCA Scatterplot") +
                  geom_text_repel(aes(label=ind_id_allIV), max.overlaps = 40)
                  
                  
                  
