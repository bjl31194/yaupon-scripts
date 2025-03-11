## 1. PACKAGES AND FILES ##

library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggrepel) 
library(readxl)
library(ggpubr)
library(rstatix)
library(patchwork)
setwd("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics")
neg_metabolites <- read_excel("./neg_metabolites.xlsx", sheet = "RT")
pos_metabolites <- read_excel("./pos_metabolites.xlsx", sheet = "RT")
quant_metabolites <- read_excel("./quant_metabolites.xlsx", sheet = "cols")


## 2. FIND DUPLICATES ##
#negative

count_na <- function(x) sum(is.na(x))  

neg_metabolites$means <- rowMeans(neg_metabolites[, 4:58], na.rm=T) 
neg_metabolites$count_na <- apply(neg_metabolites[, 4:58], 1, count_na)
neg_metabolites$is_dupe <- FALSE

neg_metabolites <- neg_metabolites %>% 
  group_by(Compound, Transition, RT)

neg_metabolites$Transition <- as.numeric(neg_metabolites$Transition)
neg_metabolites$RT <- as.numeric(neg_metabolites$RT)

# determines whether a pair of metabolites contains a duplicate using RT, Transition, and Sample Mean
for(x in 1:814) {
  if(neg_metabolites[x, 1] == neg_metabolites[x+1, 1] &
     abs(neg_metabolites[x, 2] - neg_metabolites[x+1, 2]) < 0.1 &
     abs(neg_metabolites[x, 3] - neg_metabolites[x+1, 3]) < 0.2 &
     #t.test(neg_metabolites[x, 4:55],neg_metabolites[x+1, 4:55])$p.value > 0.05) 
     abs(neg_metabolites[x, 59] - neg_metabolites[x+1, 59]) <= neg_metabolites[x, 59] * 0.10) { 
    print(neg_metabolites[x+1, 1])
    neg_metabolites[x+1, 61] <- TRUE
  }
}
#count duplicates
neg_metabolites %>%
  group_by(is_dupe) %>% tally()
#filter em out
neg_metabolites <- neg_metabolites %>%
  filter(is_dupe == FALSE) %>%
  ungroup()

newNames <- ave(as.character(neg_metabolites$Compound), 
                neg_metabolites$Compound, FUN=function(x) 
                  if (length(x)>1) paste0(x[1], '(', seq_along(x), ')') 
                else x[1]) #number compounds with the same name

neg_metabolites <- neg_metabolites %>% 
  add_column(newNames, .before = "Transition")
#get rid of unneeded cols
neg_metabolites$Compound <- NULL
colnames(neg_metabolites)[c(1)] <- c("Compound") 
neg_metabolites$Transition <- NULL
neg_metabolites$RT <- NULL
neg_metabolites$means <- NULL
neg_metabolites$count_na <- NULL
neg_metabolites$is_dupe <- NULL

# positive mode data
count_na <- function(x) sum(is.na(x))  

pos_metabolites$means <- rowMeans(pos_metabolites[, 4:58], na.rm=T) 
pos_metabolites$count_na <- apply(pos_metabolites[, 4:58], 1, count_na)
pos_metabolites$is_dupe <- FALSE

pos_metabolites <- pos_metabolites %>% 
  group_by(Compound, Transition, RT)

# determines whether a pair of metabolites contains a duplicate using RT, Transition, and Sample Mean
for(x in 1:1351) {
  if(pos_metabolites[x, 1] == pos_metabolites[x+1, 1] &
     abs(pos_metabolites[x, 2] - pos_metabolites[x+1, 2]) < 0.1 &
     abs(pos_metabolites[x, 3] - pos_metabolites[x+1, 3]) < 0.2 &
     abs(pos_metabolites[x, 59] - pos_metabolites[x+1, 59]) <= pos_metabolites[x, 59] * 0.10) { 
    print(pos_metabolites[x+1, 1])
    pos_metabolites[x+1, 61] <- TRUE
  }
}

#count duplicates
pos_metabolites %>%
  group_by(is_dupe) %>% tally()
#filter em out  
pos_metabolites <- pos_metabolites %>%
  filter(is_dupe == FALSE) %>%
  ungroup()

newNames <- ave(as.character(pos_metabolites$Compound), 
                pos_metabolites$Compound, FUN=function(x) 
                  if (length(x)>1) paste0(x[1], '(', seq_along(x), ')') 
                else x[1]) #number compounds with the same name

pos_metabolites <- pos_metabolites %>% 
  add_column(newNames, .before = "Transition")

#get rid of unneeded cols
pos_metabolites$Compound <- NULL
colnames(pos_metabolites)[c(1)] <- c("Compound") 
pos_metabolites$Transition <- NULL
pos_metabolites$RT <- NULL
pos_metabolites$means <- NULL
pos_metabolites$count_na <- NULL
pos_metabolites$is_dupe <- NULL


## 2. CHECK BLANKS ##
# negative

neg_metabolites[is.na(neg_metabolites)] <- 0.0001 #set NAs to arbitrary small value to calculate mean
blank_neg <- c()

neg_metabolites$SampleAverage <- rowMeans(neg_metabolites[, -c(1,56)])

#flags chemicals if sample average does not exceed 3 times that of the blank
for(i in 1:nrow(neg_metabolites)) {
  if(neg_metabolites[i, 57] <= 3 * neg_metabolites[i, 56]) {
    print(neg_metabolites[i, 1])
    blank_neg <- c(blank_neg, neg_metabolites[i, 1])
  }
}
cat(length(blank_neg), "chemicals found in neg blank")

#filters out blanks
neg_metabolites_filtered <- neg_metabolites[!(neg_metabolites$Compound %in% blank_neg) ,]

#cleaning up unnecessary columns
neg_metabolites_filtered[neg_metabolites_filtered == 0.0001] <- NA
neg_metabolites_filtered$SampleAverage <- NULL
neg_metabolites_filtered$blank <- NULL

write.csv(neg_metabolites_filtered, "neg_metabolites_filtered_cols.csv", row.names = FALSE)

# Check blanks - positive

pos_metabolites[is.na(pos_metabolites)] <- 0.0001 #set NAs to arbitrary small value to calculate mean
blank_pos <- c()

pos_metabolites$SampleAverage <- rowMeans(pos_metabolites[, -c(1,56)]) #get sample means

#flag if sample mean does not exceed 3 times blank level
for(i in 1:nrow(pos_metabolites)) {
  if(pos_metabolites[i, 57] <= 3 * pos_metabolites[i, 56]) {
    print(pos_metabolites[i, 1])
    blank_pos <- c(blank_pos, pos_metabolites[i, 1])
  }
}
cat(length(blank_pos), "chemicals found in pos blank")

#filter out
pos_metabolites_filtered <- pos_metabolites[!(pos_metabolites$Compound %in% blank_pos) ,]

#restore NAs
pos_metabolites_filtered[pos_metabolites_filtered == 0.0001] <- NA

#clean unnecessary cols
pos_metabolites_filtered$SampleAverage <- NULL
pos_metabolites_filtered$blank <- NULL

write.csv(pos_metabolites_filtered, "pos_metabolites_filtered_cols.csv", row.names = FALSE)


## 3. CREATE CLASS LABELS FOR SAMPLES ##
#negative
age_labels <- substr(colnames(neg_metabolites_filtered), start=4, stop=4)

treatment_labels <- substr(colnames(neg_metabolites_filtered), start=8, stop=8)

gt_labels <- substr(colnames(neg_metabolites_filtered), start=1, stop=2)

neg_metabolites_filtered_agelabels <- rbind(age_labels, neg_metabolites_filtered)
neg_metabolites_filtered_treatmentlabels <- rbind(treatment_labels, neg_metabolites_filtered)
neg_metabolites_filtered_gtlabels <- rbind(gt_labels, neg_metabolites_filtered)

neg_metabolites_filtered_agelabels[1,1] <- "Leaf Stage"
neg_metabolites_filtered_treatmentlabels[1,1] <- "Treatment"
neg_metabolites_filtered_gtlabels[1,1] <- "Genotype"

neg_metabolites_filtered_agelabels[neg_metabolites_filtered_agelabels=="M"]<-"Mature"
neg_metabolites_filtered_agelabels[neg_metabolites_filtered_agelabels=="S"]<-"Softwood"
neg_metabolites_filtered_agelabels[neg_metabolites_filtered_agelabels=="Y"]<-"Young"
neg_metabolites_filtered_treatmentlabels[neg_metabolites_filtered_treatmentlabels=="G"]<-"Green"
neg_metabolites_filtered_treatmentlabels[neg_metabolites_filtered_treatmentlabels=="R"]<-"Roasted"

write.csv(neg_metabolites_filtered_agelabels, "neg_metabolites_filtered_agelabels.csv", row.names = FALSE)
write.csv(neg_metabolites_filtered_treatmentlabels, "neg_metabolites_filtered_treatmentlabels.csv", row.names = FALSE)
write.csv(neg_metabolites_filtered_gtlabels, "neg_metabolites_filtered_gtlabels.csv", row.names = FALSE)

#make table with all metadata
neg_metabolites_filtered_alllabels <- rbind(age_labels, treatment_labels, gt_labels, neg_metabolites_filtered)
neg_metabolites_filtered_alllabels[1,1] <- "Leaf Stage"
neg_metabolites_filtered_alllabels[2,1] <- "Treatment"
neg_metabolites_filtered_alllabels[3,1] <- "Genotype"
neg_metabolites_filtered_alllabels[neg_metabolites_filtered_alllabels=="M"]<-"Mature"
neg_metabolites_filtered_alllabels[neg_metabolites_filtered_alllabels=="S"]<-"Softwood"
neg_metabolites_filtered_alllabels[neg_metabolites_filtered_alllabels=="Y"]<-"Young"
neg_metabolites_filtered_alllabels[neg_metabolites_filtered_alllabels=="G"]<-"Green"
neg_metabolites_filtered_alllabels[neg_metabolites_filtered_alllabels=="R"]<-"Roasted"

#positive

age_labels <- substr(colnames(pos_metabolites_filtered), start=4, stop=4)

treatment_labels <- substr(colnames(pos_metabolites_filtered), start=8, stop=8)

gt_labels <- substr(colnames(pos_metabolites_filtered), start=1, stop=2)

pos_metabolites_filtered_agelabels <- rbind(age_labels, pos_metabolites_filtered)
pos_metabolites_filtered_treatmentlabels <- rbind(treatment_labels, pos_metabolites_filtered)
pos_metabolites_filtered_gtlabels <- rbind(gt_labels, pos_metabolites_filtered)

pos_metabolites_filtered_agelabels[1,1] <- "Leaf Stage"
pos_metabolites_filtered_treatmentlabels[1,1] <- "Treatment"
pos_metabolites_filtered_gtlabels[1,1] <- "Genotype"

pos_metabolites_filtered_agelabels[pos_metabolites_filtered_agelabels=="M"]<-"Mature"
pos_metabolites_filtered_agelabels[pos_metabolites_filtered_agelabels=="S"]<-"Softwood"
pos_metabolites_filtered_agelabels[pos_metabolites_filtered_agelabels=="Y"]<-"Young"
pos_metabolites_filtered_treatmentlabels[pos_metabolites_filtered_treatmentlabels=="G"]<-"Green"
pos_metabolites_filtered_treatmentlabels[pos_metabolites_filtered_treatmentlabels=="R"]<-"Roasted"

write.csv(pos_metabolites_filtered_agelabels, "pos_metabolites_filtered_agelabels.csv", row.names = FALSE)
write.csv(pos_metabolites_filtered_treatmentlabels, "pos_metabolites_filtered_treatmentlabels.csv", row.names = FALSE)
write.csv(pos_metabolites_filtered_gtlabels, "pos_metabolites_filtered_gtlabels.csv", row.names = FALSE)

# create class labels - quant

age_labels <- substr(colnames(quant_metabolites), start=4, stop=4)

treatment_labels <- substr(colnames(quant_metabolites), start=8, stop=8)

quant_metabolites_agelabels <- rbind(age_labels, quant_metabolites)
quant_metabolites_treatmentlabels <- rbind(treatment_labels, quant_metabolites)


quant_metabolites_agelabels[1,1] <- "Age"
quant_metabolites_treatmentlabels[1,1] <- "Treatment"

write.csv(quant_metabolites_agelabels, "quant_metabolites_agelabels.csv", row.names = FALSE)
write.csv(quant_metabolites_treatmentlabels, "quant_metabolites_treatmentlabels.csv", row.names = FALSE)


## 4. POST-METABOANALYST: FIND COMPOUNDS DIFFERENTIATED BY LEAF STAGE ##
## using overlap b/t linear model/ANOVA and random forest results

# negative mode data
neg_lm_results_age <- read.csv("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics/Neg metabolite analysis/multivariate//lm_age_all.csv")

neg_rf_sigfeatures_age <- read.csv("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics/Neg metabolite analysis/multivariate//rf_sigfeatures_age.csv")

# Make table with both lm results and rf
lm_sorted <- neg_lm_results_age %>%
  arrange(X)
rf_sorted <- neg_rf_sigfeatures_age %>%
  arrange(X)
identical(lm_sorted$X, rf_sorted$X)
lm_sorted$MeanDecreaseAccuracy <- rf_sorted$MeanDecreaseAccuracy 
neg_lm_sorted <- lm_sorted %>%
  arrange(adj.P.Val)
colnames(neg_lm_sorted)[1] <- "Compound" 
colnames(neg_lm_sorted)[8] <- "MeanDecreaseAccuracy"

#filter by FC and adjusted p value
neg_lm_sigup <- neg_lm_sorted %>% 
  filter(mature.young > 2, adj.P.Val < 0.05, MeanDecreaseAccuracy > quantile(rf_sigfeatures_age$MeanDecreaseAccuracy, 0.75)) 

neg_lm_sigdown <- neg_lm_sorted %>% 
  filter(mature.young < -2, adj.P.Val < 0.05, MeanDecreaseAccuracy > quantile(rf_sigfeatures_age$MeanDecreaseAccuracy, 0.75))

write.csv(neg_lm_sigup, "/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics/Neg metabolite analysis/multivariate/neg_sigup_with_age.csv", row.names = FALSE)
write.csv(neg_lm_sigdown, "/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics/Neg metabolite analysis/multivariate/neg_sigdown_with_age.csv", row.names = FALSE)

## Positive mode data 

setwd("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics/Pos metabolite analysis/multivariate/")

pos_lm_results_age <- read.csv("lm_age_all.csv")

pos_rf_sigfeatures_age <- read.csv("rf_sigfeatures_age.csv")

# Make table with both lm results and rf
lm_sorted <- pos_lm_results_age %>%
  arrange(X)
rf_sorted <- pos_rf_sigfeatures_age %>%
  arrange(X)
identical(lm_sorted$X, rf_sorted$X)
lm_sorted$MeanDecreaseAccuracy <- rf_sorted$MeanDecreaseAccuracy 
pos_lm_sorted <- lm_sorted %>%
  arrange(adj.P.Val)
colnames(pos_lm_sorted)[1] <- "Compound" 
colnames(pos_lm_sorted)[8] <- "MeanDecreaseAccuracy"

#filter by FC and adjusted p value
pos_lm_sigup <- pos_lm_sorted %>% 
  filter(mature.young > 2, adj.P.Val < 0.05, MeanDecreaseAccuracy > quantile(rf_sigfeatures_age$MeanDecreaseAccuracy, 0.75)) 

pos_lm_sigdown <- pos_lm_sorted %>% 
  filter(mature.young < -2, adj.P.Val < 0.05, MeanDecreaseAccuracy > quantile(rf_sigfeatures_age$MeanDecreaseAccuracy, 0.75))

write.csv(pos_lm_sigup, "/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics/Pos metabolite analysis/multivariate/pos_sigup_with_age.csv", row.names = FALSE)
write.csv(pos_lm_sigdown, "/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics/Pos metabolite analysis/multivariate/pos_sigdown_with_age.csv", row.names = FALSE)

#summary
cat(nrow(neg_sigup_with_age_overlap), "chemicals up in neg mode")
cat(nrow(neg_sigdown_with_age_overlap), "chemicals down in neg mode")
cat(nrow(pos_sigup_with_age_overlap), "chemicals up in pos mode")
cat(nrow(pos_sigdown_with_age_overlap), "chemicals down in pos mode")


## 5. DATA VIS / FIGURE GENERATION ##

## Targeted Metabolomics ##

# create class labels for samples
quant_metabolites <- read_excel("./quant_metabolites.xlsx", sheet = "rows")

age_labels <- substr(quant_metabolites$Sample, start=4, stop=4)

treatment_labels <- substr(quant_metabolites$Sample, start=8, stop=8)

gt_labels <- substr(quant_metabolites$Sample, start=1, stop=2)

quant_metabolites <- bind_cols(age_labels, quant_metabolites)
quant_metabolites <- bind_cols(treatment_labels, quant_metabolites)
quant_metabolites <- bind_cols(gt_labels, quant_metabolites)


quant_metabolites[quant_metabolites=="M"] <- "mature"
quant_metabolites[quant_metabolites=="S"] <- "softwood"
quant_metabolites[quant_metabolites=="Y"] <- "young"
quant_metabolites[quant_metabolites=="G"] <- "green"
quant_metabolites[quant_metabolites=="R"] <- "roasted"

colnames(quant_metabolites)[2] <- "Treatment"
colnames(quant_metabolites)[3] <- "Age"
colnames(quant_metabolites)[1] <- "Genotype"

quant_metabolites$Age <- factor(quant_metabolites$Age, levels = c("young", "softwood", "mature"))

# get group stats
aggregate(quant_metabolites[, 5:8], list(quant_metabolites$Age), range)

## t tests for treatment and genotype effect ##

#caffeine test
group_by(quant_metabolites, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(Caffeine, na.rm = TRUE),
    sd = sd(Caffeine, na.rm = TRUE)
  )
var.test(Caffeine ~ Treatment, data = quant_metabolites)

caf <- t.test(Caffeine ~ Treatment, data = quant_metabolites, var.equal = TRUE)
caf

#theobromine test
group_by(quant_metabolites, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(Theobromine, na.rm = TRUE),
    sd = sd(Theobromine, na.rm = TRUE)
  )
var.test(Theobromine ~ Treatment, data = quant_metabolites)

tb <- t.test(Theobromine ~ Treatment, data = quant_metabolites, var.equal = TRUE)
tb

#theacrine test
group_by(quant_metabolites, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(Theacrine, na.rm = TRUE),
    sd = sd(Theacrine, na.rm = TRUE)
  )
var.test(Theacrine ~ Treatment, data = quant_metabolites)

caf <- t.test(Theacrine ~ Treatment, data = quant_metabolites, var.equal = TRUE)
caf

#CGA test
group_by(quant_metabolites, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(CGA, na.rm = TRUE),
    sd = sd(CGA, na.rm = TRUE)
  )
var.test(CGA ~ Treatment, data = quant_metabolites)

caf <- t.test(CGA ~ Treatment, data = quant_metabolites, var.equal = TRUE)
caf

#caffeine test - GT
group_by(quant_metabolites, Genotype) %>%
  summarise(
    count = n(),
    mean = mean(Caffeine, na.rm = TRUE),
    sd = sd(Caffeine, na.rm = TRUE)
  )

caf <- aov(Caffeine ~ Genotype, data = quant_metabolites) 
summary(caf)
tb <- aov(Theobromine ~ Genotype, data = quant_metabolites) 
summary(tb)
tc <- aov(Theacrine ~ Genotype, data = quant_metabolites) 
summary(tc)
cga <- aov(CGA ~ Genotype, data = quant_metabolites) 
summary(cga)
TukeyHSD(tc)

#theobromine test - GT
group_by(quant_metabolites, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(Theobromine, na.rm = TRUE),
    sd = sd(Theobromine, na.rm = TRUE)
  )
var.test(Theobromine ~ Treatment, data = quant_metabolites)

tb <- t.test(Theobromine ~ Treatment, data = quant_metabolites, var.equal = TRUE)
tb

#theacrine test - GT
group_by(quant_metabolites, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(Theacrine, na.rm = TRUE),
    sd = sd(Theacrine, na.rm = TRUE)
  )
var.test(Theacrine ~ Treatment, data = quant_metabolites)

caf <- t.test(Theacrine ~ Treatment, data = quant_metabolites, var.equal = TRUE)
caf

#CGA test - GT
group_by(quant_metabolites, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(CGA, na.rm = TRUE),
    sd = sd(CGA, na.rm = TRUE)
  )
var.test(CGA ~ Treatment, data = quant_metabolites)

caf <- t.test(CGA ~ Treatment, data = quant_metabolites, var.equal = TRUE)
caf

## Targeted metabolomics by age grouped by treatment ##
my_comparisons <- list(c("mature", "softwood"), c("softwood", "young"), c("mature", "young"))

my_theme <- theme(
  axis.title.x = element_text(size = 11, face="bold"),
  axis.text.x = element_text(size = 9, face="bold"),
  axis.title.y = element_text(size = 11, face="bold"),
  axis.text.y = element_text(size = 9, face="bold"),
  legend.text = element_text(size = 10),
  legend.title = element_text(size=11),
  plot.title = element_text(size = 14, face = "bold"))
fs <- 4

#individual plots
caffeine <- quant_metabolites %>%
  ggplot(aes(x=Age, y=Caffeine, fill=Treatment)) + 
  geom_boxplot() +
  ggtitle("Caffeine") +
  xlab(" ") + ylab("mg/g") +
  scale_fill_manual(values = c("forestgreen","tan3")) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = .1)) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
caffeine$layers[[2]]$aes_params$textsize <- fs
#caffeine

theobromine <- quant_metabolites %>%
  ggplot(aes(x=Age, y=Theobromine, fill=Treatment)) + 
  geom_boxplot() +
  ggtitle("Theobromine") +
  xlab(" ") + ylab(" ") +
  scale_fill_manual(values = c("forestgreen","tan3")) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = .1)) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
theobromine$layers[[2]]$aes_params$textsize <- fs
#theobromine

theacrine <- quant_metabolites %>%
  ggplot(aes(x=Age, y=Theacrine, fill=Treatment)) + 
  geom_boxplot() +
  ggtitle("Theacrine") +
  xlab("Leaf Stage") + ylab("mg/g") +
  scale_fill_manual(values = c("forestgreen","tan3")) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = .1)) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
theacrine$layers[[2]]$aes_params$textsize <- fs
#theacrine

cga <- quant_metabolites %>%
  ggplot(aes(x=Age, y=CGA, fill=Treatment)) + 
  geom_boxplot() +
  ggtitle("Chlorogenic Acid") +
  xlab("Leaf Stage") + ylab(" ") +
  scale_fill_manual(values = c("forestgreen","tan3")) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = .1)) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
cga$layers[[2]]$aes_params$textsize <- fs
#cga

(caffeine | theobromine) / #grouped plot
(theacrine | cga) +
  plot_layout(guides = "collect", axes = "collect_x") + 
  plot_annotation(tag_levels = "A") +
    theme(legend.position = 'bottom')
  
## Targeted metabolomics by genotype ##

AR <- quant_metabolites[quant_metabolites$Sample %like% "AR" ,]
VA <- quant_metabolites[quant_metabolites$Sample %like% "VA" ,]
FL <- quant_metabolites[quant_metabolites$Sample %like% "FL" ,]

#set this variable to select which genotype to plot from
gt <- FL

#individual plots - 
caffeineGT <- gt %>%
  ggplot(aes(x=Age, y=Caffeine, fill=Treatment)) + 
  geom_boxplot() +
  ggtitle("Caffeine") +
  xlab(" ") + ylab("mg/g") +
  scale_fill_manual(values = c("forestgreen","tan3")) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = .1)) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
caffeineGT$layers[[2]]$aes_params$textsize <- fs
caffeineGT

theobromineGT <- gt %>%
  ggplot(aes(x=Age, y=Theobromine, fill=Treatment)) + 
  geom_boxplot() +
  ggtitle("Theobromine") +
  xlab(" ") + ylab(" ") +
  scale_fill_manual(values = c("forestgreen","tan3")) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = .1)) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
theobromineGT$layers[[2]]$aes_params$textsize <- fs
theobromineGT

theacrineGT <- gt %>%
  ggplot(aes(x=Age, y=Theacrine, fill=Treatment)) + 
  geom_boxplot() +
  ggtitle("Theacrine") +
  xlab("Leaf Stage") + ylab("mg/g") +
  scale_fill_manual(values = c("forestgreen","tan3")) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = .1)) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
theacrineGT$layers[[2]]$aes_params$textsize <- fs
theacrineGT

cgaGT <- gt %>%
  ggplot(aes(x=Age, y=CGA, fill=Treatment)) + 
  geom_boxplot() +
  ggtitle("Chlorogenic Acid") +
  xlab("Leaf Stage") + ylab(" ") +
  scale_fill_manual(values = c("forestgreen","tan3")) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = .1)) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
cgaGT$layers[[2]]$aes_params$textsize <- fs
cgaGT

(caffeineGT | theobromineGT) / #make grouped plots
  (theacrineGT | cgaGT) +
  plot_layout(guides = "collect", axes = "collect_x") + 
  plot_annotation(tag_levels = "A") +
  theme(legend.position = 'bottom')

## Volcano plots - NEG - TREATMENT ##
setwd("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics/Neg metabolite analysis/univariate-by-treatment/")
neg <- read_excel("neg_volcano_all.xlsx")

#clean up
colnames(neg)[c(1,3,4,5)] <- c("Compound", "log2FC", "p.adj", "-log10(pval)")
neg$Compound[neg$Compound == "1-O-Feruloyl-Œ≤-D-glucose"] <- "1-O-Feruloyl-beta-D-glucose"

#create column for sig compounds
neg$diff[neg$log2FC > 1 & neg$p.adj < 0.05] <- "UP"
neg$diff[neg$log2FC < -1 & neg$p.adj < 0.05] <- "DOWN"

write_excel_csv(neg, "neg_sig_treatment.csv")

#select labeled compounds - speciic list or top X
neg$label <- ifelse(neg$Compound %in% c("Asiaticoside(1)","4-O-Caffeoylshikimic acid(3)","4-O-Caffeoylshikimic acid(1)", "Puerarin xyloside(1)","5-O-Feruloylquinic acid(2)","Soyasaponin I(1)","Pseudopurpurin","6''-O-Acetylgenistin","Enicoflavine","Larixinic Acid","1-O-Feruloyl-beta-D-glucose"), neg$Compound, NA)
#neg$label <- ifelse(neg$Compound %in% head(neg[order(neg$p.adj), "Compound"], 10), neg$Compound, NA)

vol_theme = theme( #this theme will be used for both volcano plots
  axis.title.x = element_text(size = 14, face="bold"),
  axis.text.x = element_text(size = 12, face="bold"),
  axis.title.y = element_text(size = 14, face="bold"),
  axis.text.y = element_text(size = 12, face="bold"),
  legend.text = element_text(size = 12),
  legend.title = element_text(size=12),
  plot.title = element_text(size = 18, face = "bold"))
fs <- 4

volcano_neg <- neg %>% #create plot
  ggplot(aes(x = log2FC, y = `-log10(pval)`, col = diff, label = label)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("forestgreen", "tan3", "grey"), # to set the colours of our variable
                     labels = c("Down in roasted", "Up in roasted", "Not significant")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 12), xlim = c(-3, 7)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = '', #legend_title
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  scale_x_continuous(breaks = seq(-3, 7, 1)) + # to customise the breaks in the x axis
  ggtitle('Chemicals in negative mode significantly affected by the roasting process') + # Plot title
  geom_text_repel(size = 6, max.overlaps = Inf, show.legend = FALSE) + # To show all labels 
  vol_theme +
  theme(legend.position="none")
volcano_neg  

## Volcano plots - POS - TREATMENT ##

setwd("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics/Pos metabolite analysis/univariate-by-treatment/")
pos <- read_excel("pos_volcano_all.xlsx")

#set column names and clean up
colnames(pos)[c(1,3,4,5)] <- c("Compound", "log2FC", "p.adj", "-log10(pval)") 
pos$Compound[pos$Compound == "1,1'-(Tetrahydro-6a-hydroxy-2,3a,5-trimethylfuro[2,3-d]-1,3-dioxole-2,5-diyl)bis-ethanone"] <- "Diacetyl trimer"

#create column for sig compounds
pos$diff[pos$log2FC > 1 & pos$p.adj < 0.05] <- "UP"
pos$diff[pos$log2FC < -1 & pos$p.adj < 0.05] <- "DOWN"

write_excel_csv(pos, "pos_sig_treatment.csv")

#pos$label <- ifelse(pos$Compound %in% head(pos[order(pos$p.adj), "Compound"], 10), pos$Compound, NA)
pos$label <- ifelse(pos$Compound %in% c("Isonicotineamide","Ecgonine","N-butanoyl-lhomoserine lactone","Pentanamide","Quinone(2)","Diacetyl trimer","5-O-Caffeoylshikimic acid(7)","Quercetin(3)","D-1-[(3-Carboxypropyl)amino]-1-deoxyfructose","5-O-Feruloylquinic acid(3)","Cepharadione B","3,4-Dicaffeoyl-1,5-quinolactone(4)","Lappaol A","3-Demethylsimmondsin 2'-(Z)-ferulate","5,7alpha-Dihydro-1,4,4,7a-tetramethyl-4H-indene"), pos$Compound, NA)

volcano_pos <- pos %>% #make volcano plot
  ggplot(aes(x = log2FC, y = `-log10(pval)`, col = diff, label = label)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("forestgreen", "tan3", "grey"), # to set variable colors 
                     labels = c("Down in roasted", "Up in roasted", "Not significant")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 20), xlim = c(-4, 8)) + 
  labs(color = '', #legend_title
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  scale_x_continuous(breaks = seq(-3, 7, 1)) + # to customise the breaks in the x axis
  ggtitle('Chemicals in positive mode significantly affected by the roasting process') + # Plot title
  geom_text_repel(size = 6, max.overlaps = Inf, show.legend = FALSE) + # To show all labels 
  vol_theme
volcano_pos  

#plot two together
(volcano_neg | volcano_pos) + 
  plot_layout(guides = "collect", axes = "collect_y") + 
  plot_annotation(tag_levels = "A")

#find overlaps in lists
up_pos <- pos %>% 
  filter(grepl('UP', diff))
    
up_neg <- neg %>%
  filter(grepl('UP', diff))

up <- inner_join(up_pos, up_neg, by = c('Compound'))


down_pos <- pos %>% 
  filter(grepl('DOWN', diff))

down_neg <- neg %>%
  filter(grepl('DOWN', diff))

down <- inner_join(down_pos, down_neg, by = c('Compound'))

allcomps <- inner_join(pos, neg, by = c('Compound'))
up_both <- inner_join(neg_sigup_with_age_overlap, pos_sigup_with_age_overlap, by = c('X'))
down_both <- inner_join(neg_sigdown_with_age_overlap, pos_sigdown_with_age_overlap, by = c('X'))
lm_up_both <- inner_join(neg_lm_sigup, pos_lm_sigup, by = c('X'))
lm_down_both <- inner_join(neg_lm_sigdown, pos_lm_sigdown, by = c('X'))
