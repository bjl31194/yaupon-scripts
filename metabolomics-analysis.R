#libraries & data import

library(data.table)
library(tidyverse)
library(ggplot2)
setwd("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics")
neg_metabolites <- read.csv("./neg-metabolites.csv")

# Check blank

neg_metabolites[is.na(neg_metabolites)] <- 0.0001
blank <- c()

for(i in 2:length(neg_metabolites)) {
  sampleAverage <- mean(neg_metabolites[1:54, i])
  if(sampleAverage <= 3 * neg_metabolites[55, i]) {
    print(colnames(neg_metabolites[i]))
    blank <- c(blank, colnames(neg_metabolites[i]))
  }
}

neg_metabolites_filtered <- neg_metabolites[ , !(names(neg_metabolites) %in% blank)]

neg_metabolites_filtered[neg_metabolites_filtered == 0.0001] <- NA

write.csv(neg_metabolites_filtered, "neg_metabolites_filtered.csv", row.names = FALSE)

#Normalization using average internal standard peak area ratio

rownames(neg_metabolites) <- neg_metabolites$Sample

neg_metabolites$Sample <- NULL

slice = neg_metabolites[, 6:818]

IS <- t(neg_metabolites[5])
IS[55] <- 1

norm <- slice * IS

norm <- cbind(neg_metabolites[, c(1,3)], norm)

norm[is.na(norm)] <- 0.0001

present_in_blank <- row.names(as.data.frame(grep(0.0001, norm, value = TRUE, invert = TRUE)))

write_csv(present_in_blank, "neg_present_in_blank.csv")
## sort rows by roasted or green, plus age ##

greenrows <- grep("-G", row.names(norm))

greenDF <- norm[greenrows ,]

normSorted <- norm[greenrows ,]

roastedrows <-grep("-R", row.names(norm))

roastedDF <- norm[roastedrows ,]

normSorted <- rbind(roastedDF, normSorted)

normSorted[, 'Treatment'] <- NA

normSorted <- normSorted %>% relocate(Treatment, .before = colnames(normSorted[1]))

normSorted[1:27,1] <- "roasted"
normSorted[28:54,1] <- "green"

youngrows <- normSorted[grep("-Y-", row.names(normSorted)) ,]

softrows <- normSorted[grep("-S-", row.names(normSorted)) ,]

matrows <- normSorted[grep("-M-", row.names(normSorted)) ,]

normSorted <- rbind(youngrows, softrows, matrows)

normSorted[, 'Age'] <- NA

normSorted <- normSorted %>% relocate(Age, .before = colnames(normSorted[1]))

normSorted[1:18,1] <- "young"
normSorted[19:36,1] <- "softwood"
normSorted[37:54,1] <- "mature"


# # Pareto Scaling and PCA # #
install.packages("IMIFA")
library(IMIFA)
install.packages("factoextra")
library(factoextra)

normPareto <- pareto_scale(normSorted[, 3:815], centering=TRUE)
pca <- prcomp(normPareto, scale = FALSE)
fviz_eig(pca)
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             select.var = list(contrib=20)
)

fviz_contrib(pca,
             choice = "var",
             axes = 2,
             top = 10, color = 'darkorange3', barfill  = 'blue4',fill ='blue4')

fviz_pca_ind(pca,
             label = "none", # hide individual labels
             habillage = normSorted$Treatment, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE # Concentration ellipses
)

fviz_pca_biplot(pca, 
                repel = TRUE, 
                select.var = list(contrib = 10)
)

a <- pca$rotation

a %>% as.data.frame %>% rownames_to_column %>% 
  select(rowname, PC1, PC2) %>% arrange(desc(PC1^2+PC2^2)) %>% head(10)

# # t test and ANOVA # #

t_test <- lapply(normSorted[3:815], function(x) t.test(x~normSorted$Treatment)$p.value)

fdrs <- p.adjust(t_test, method="BH")

bonf <- p.adjust(t_test, method="bonferroni")

fdrs_sorted <- sort(fdrs)

bonf_sorted <- sort(bonf)

anova <- lapply(normSorted[3:815], function(x) aov(x~normSorted$Age))

#pvals_anova <- lapply(normSorted[3:815], function(x) summary(aov(x~normSorted$Age))[[1]]$`Pr(>F)`[1])

for (i in 3:815) {
  column <- names(normSorted[i])
  anova <- broom::tidy(aov(normSorted[,i] ~ Age, data = normSorted))
  
  # only want aov with P < 0.05 printed
  if(anova$p.value[1] < 0.05) {
    
    print(column)
    print(anova)
  }
}

# post-hoc tests

tukey_results <- list(NULL)

for (i in 3:815) {
  column <- names(normSorted[i])
  anova <- aov(normSorted[,i] ~ Age, data = normSorted)
  tukey <- TukeyHSD(anova, conf.level=0.99)
  
  # only want tukey with P < 0.05 printed
  if(any(tukey$Age[, "p adj"] < 0.01)) {
    j <- i-2
    tukey_results[j] <- tukey
  }
}


# only want tukey with P < 0.05 printed
if(any(tukey$treatment[, "p adj"] < 0.05)) {
  
  print(column)
  print(setNames(tukey, column))
}


#get average peak areas for all chemicals in roasted and green

greenMeans <- sapply(greenDF, mean)

roastedMeans <- sapply(roastedDF, mean)

# get ratios of chemical levels in green vs. roasted leaves

RG_ratios <- data.frame(greenMeans, roastedMeans, roastedMeans/greenMeans)
RG_ratios <- data.frame(row.names(RG_ratios), greenMeans, roastedMeans, roastedMeans/greenMeans)
colnames(RG_ratios)[4] = "Ratio"
colnames(RG_ratios)[1] = "Chemical"

RG_ratios_sorted <- RG_ratios %>%
  arrange(desc(`Ratio`))

# plot peak areas of top ten chemicals enriched in green and roasted treatments

topRoasted <- head(RG_ratios_sorted, 10)

plotRoasted <- pivot_longer(topRoasted, cols = c(greenMeans,roastedMeans), names_to = "Treatment", values_to = "Peak_Area")

topRoastedPlot <- ggplot(plotRoasted, aes(x = Chemical, y = Peak_Area, fill = Treatment)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("greenMeans" = "forestgreen",
                               "roastedMeans" = "tan3")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_flip() 
 
  
topRoastedPlot

topGreen <- tail(RG_ratios_sorted, 10)

plotGreen <- pivot_longer(topGreen, cols = c(greenMeans,roastedMeans), names_to = "Treatment", values_to = "Peak_Area")

topGreenPlot <- ggplot(plotGreen, aes(x = Chemical, y = Peak_Area, fill = Treatment)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("greenMeans" = "forestgreen",
                               "roastedMeans" = "tan3")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_flip() 
topGreenPlot

# plot R/G log ratios

top20 <- rbind(topRoasted, topGreen)

top20$logRatio <- log(top20$Ratio)

plot20 <- pivot_longer(top20, cols = c(greenMeans,roastedMeans), names_to = "Treatment", values_to = "Peak_Area")


top20Plot <- top20 %>%
  mutate(Chemical = fct_reorder(Chemical, logRatio)) %>%
  ggplot( aes(x = Chemical, y = logRatio)) +
  geom_point(aes(color = logRatio > 0),
             show.legend = FALSE) +
  scale_color_manual(values = c("forestgreen","tan3")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed")
top20Plot

  ## excluding Enicoflavine 
top20Plot <- top20 %>%
  mutate(Chemical = fct_reorder(Chemical, logRatio)) %>%
  filter(logRatio < 20) %>%
  ggplot( aes(x = Chemical, y = logRatio)) +
  geom_point(aes(color = logRatio > 0),
             show.legend = FALSE) +
  scale_color_manual(values = c("forestgreen","tan3")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed")
top20Plot

## sort based on tissue age

youngrows <- norm[grep("-Y-", row.names(norm)) ,]

softrows <- norm[grep("-S-", row.names(norm)) ,]

matrows <- norm[grep("-M-", row.names(norm)) ,]

normSortedAge <- rbind(youngrows, softrows, matrows)

#normSortedAgeLog <- log(normSortedAge)

## calculate means 

youngMeans <- sapply(youngrows, mean)
softMeans <- sapply(softrows, mean)
matMeans <- sapply(matrows, mean)

Age <- data.frame(youngMeans, softMeans, matMeans)
Age <- data.frame(row.names(Age), youngMeans, softMeans, matMeans)
colnames(Age)[1] = "Chemical"

## create qualifying variable to classify based on tissue age

normSortedAge[, 'Age'] <- NA

normSortedAge <- normSortedAge %>% relocate(Age, .before = colnames(normSortedAge[1]))

normSortedAge[1:18,1] <- "young"
normSortedAge[19:36,1] <- "softwood"
normSortedAge[37:54,1] <- "mature"

## ANOVA for tissue types

test <- normSortedAge %>% 
  select(Age, Malonic.acid_105.0205)

summary(test)
ggplot(test) +
  aes(x = Age, y = Malonic.acid_105.0205, color = Age) +
  geom_jitter() +
  theme(legend.position = "none")

res_aov <- aov(Malonic.acid_105.0205 ~ Age,
               data = test
            )
par(mfrow = c(1, 2)) # combine plots

## histogram
hist(res_aov$residuals)

## QQ-plot
library(car)
qqPlot(res_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)

leveneTest(Malonic.acid_105.0205 ~ Age,
           data = test
)

summary(res_aov)

# Trim molecular weights from list

chemicals <- read.csv("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/metabolomics/Pos metabolite analysis/univariate/by_treatment/Download (2)/volcano_pos_treatment.csv")

names <- chemicals[, 1]

names_shortened <- as.data.frame(gsub("\\_.*","",names))

write.csv(names_shortened,"./pos_volcano_treatment_names.csv", row.names = TRUE)



