################
## Packages ##
################
install.packages("tidypopgen")
install.packages("ggtext")
library(paletteer)
require(tidypopgen)
library(rnaturalearth)
library(ggplot2)
library(tidyverse)
library(ggtext)
library(LEA)
library(ggrepel)
library(sf)
setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1-5")

'%!in%' <- function(x,y)!('%in%'(x,y))
###########################
## Data import ##
###########################

# import data from VCF to gen tibble
Ilex_gt <- gen_tibble("Ilex1-5_vcftools_filter.vcf.gz")
Ilex_gt <- gt_load('/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1-5/Ilex1-5_vcftools_filter.vcf.gt')
Ivom_ig_gt <- gen_tibble("Ivom1-5_filter_names.ann.intergenic.vcf")
#Ivom_ig_gt <- gt_load("Ivom1-5_filter_names.ann.intergenic.gt")
# add population metadata
pops <- read.csv("Ilex_meta.csv")
wild_pops <- read.csv("Ivom1-5_wild_meta.csv")
pops <- as_tibble(pops)
Ilex_gt <- Ilex_gt %>% left_join(pops, by = "id")
#Ivom_ig_gt <- Ivom_ig_gt %>% left_join(pops, by = "id")
Ilex_gt %>% glimpse()

# visualize on map
Ivom_gt <- gt_add_sf(Ivom_gt, c("lon", "lat"))
map <- ne_countries(
  continent = "North America",
  type = "map_units", scale = "medium"
)

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = Ivom_gt$geometry) +
  coord_sf(
    xlim = c(-105, -75),
    ylim = c(20, 42)
  ) +
  theme_minimal()

###############################
## for use with RDA script ##
###############################

# prepare for export to RDA
# non-cultivated IVV only 
flagged_feral <- c("VA-FR-2","VA-FR-2X","FL-BE-1","VA-FR-1","VA-FR-1X","NC-CB-1","VA-FR-3","VA-FR-3X","NC-KW-2","NC-CB-2")
Ivom_wild_gt <- Ivom_gt %>%
  filter(population != "cultivated") %>%
  filter(id %!in% flagged_feral) %>%
  filter(site != "MC-AR")
Ivom_wild_gt <- Ivom_wild_gt %>% select_loci_if(loci_maf(genotypes) > 0)
#Ivom_wild_gt <- Ivom_wild_gt %>% filter(indiv_missingness(genotypes) < 0.2)
Ivom_wild_gt <- Ivom_wild_gt %>% select_loci_if(loci_missingness(genotypes) < 0.1)
Ivom_wild_gt <- Ivom_wild_gt %>% select_loci_if(loci_maf(genotypes) > 0.05)
Ivom_wild_gt <- gt_update_backingfile(Ivom_wild_gt, backingfile = "Ivom_wild_gt.rds")
Ivom_wild_gt <- gt_impute_simple(Ivom_wild_gt, method = "mode")
Ivom_wild_gt <- Ivom_wild_gt %>%
  select_loci_if(loci_ld_clump(Ivom_wild_gt, thr_r2 = 0.2))
gt_uses_imputed(Ivom_wild_gt)
gt_set_imputed(Ivom_wild_gt, TRUE)
gen <- gt_as_genind(Ivom_wild_gt)

# prepare intergenic dataset 
Ivom_ig_gt <- Ivom_ig_gt %>%
  filter(id %in% Ivom_wild_gt$id)
#Ivom_ig_gt <- Ivom_ig_gt %>% filter(indiv_missingness(genotypes) < 0.2)
Ivom_ig_gt <- Ivom_ig_gt %>% select_loci_if(loci_missingness(genotypes) < 0.1)
#Ivom_ig_gt <- Ivom_ig_gt %>% select_loci_if(loci_maf(genotypes) > 0.01)
Ivom_ig_gt <- gt_update_backingfile(Ivom_ig_gt, backingfile = tempfile())
Ivom_ig_gt <- gt_impute_simple(Ivom_ig_gt, method = "mode")
gt_uses_imputed(Ivom_ig_gt)
gt_set_imputed(Ivom_ig_gt, TRUE)

# 'gen' is used as input for RDA in different script 

###############
## Data QC ##
###############
# for full Ilex data
indiv_qc_Ilex <- Ilex_gt %>% qc_report_indiv()
autoplot(indiv_qc_Ilex, type = "scatter")
loci_qc_Ilex <- Ilex_gt %>% qc_report_loci()
autoplot(loci_qc_Ivom, type = "maf")
autoplot(loci_qc_Ivom, type = "missing")

## for Ivom data
Ivom_gt <- Ilex_gt %>%
  filter(spp == "vomitoria")
# remove monomorphic loci
Ivom_gt <- Ivom_gt %>% select_loci_if(loci_maf(genotypes) > 0)

## prepare data for gemma ##
sexed_gt <- Ilex_gt %>%
  filter(spp == "vomitoria") %>%
  filter(population != "cultivated") %>%
  filter(id %!in% flagged_feral) %>%
  filter(site != "MC-AR") %>%
  filter(sex == "M" | sex == "F")
sexed_gt <- sexed_gt %>% select_loci_if(loci_maf(genotypes) > 0.05)
sexed_gt <- sexed_gt %>% select_loci_if(loci_missingness(genotypes) < 0.1)
sexed_gt <- gt_update_backingfile(sexed_gt, backingfile = "sexed_gt.rds")
sexed_gt <- gt_impute_simple(sexed_gt, method = "mode")

gt_as_plink(sexed_gt, file = "Ivom_wild_sexed.bed", type = "bed")
sexed_gt %>%
  select(id,sex) %>%
  write_tsv("Ivom_wild_sex_phenotypes.txt")

## prepare dataset for rangeExpansion ##

rp_gt <- Ilex_gt %>%
  filter(spp == "vomitoria" | spp == "mxsubspp")

rp_gt <- rp_gt %>%
  filter(population != "cultivated") %>%
  filter(id %!in% flagged_feral) %>%
  filter(site != "MC-AR")
rp_gt <- rp_gt %>% select_loci_if(loci_maf(genotypes) > 0)
rp_gt <- rp_gt %>% select_loci_if(loci_maf(genotypes) > 0.05)
rp_gt <- rp_gt %>% select_loci_if(loci_missingness(genotypes) < 0.1)
rp_gt <- gt_update_backingfile(rp_gt, backingfile = "rp_gt.rds")
rp_gt <- gt_impute_simple(rp_gt, method = "mode")
rp_gt <- rp_gt %>%
  select_loci_if(loci_ld_clump(rp_gt, thr_r2 = 0.2))

gt_as_plink(rp_gt, file = "Ivom_wild_rp_out_all.bed", type = "bed")
rp_gt %>%
  select(id,lat,lon,site) %>%
  write_tsv("RangeExp_coords_all.tsv")
#check out data by individuals, QC
Ivom_gt %>% show_genotypes(indiv_indices = 1:5, loci_indices = 1:10)
head(Ivom_gt %>% show_loci())
indiv_qc_Ivom <- Ivom_gt %>% qc_report_indiv()
autoplot(indiv_qc_Ivom, type = "scatter")
#filter out individuals with >20% missing data
Ivom_gt <- Ivom_gt %>% filter(indiv_missingness(genotypes) < 0.2)
#Ivom_ig_gt <- Ivom_ig_gt %>% filter(indiv_missingness(genotypes) < 0.2)
# check out data by loci
Ivom_gt <- Ivom_gt %>% group_by(site)
loci_qc_Ivom <- Ivom_gt %>% qc_report_loci()
autoplot(loci_qc_Ivom, type = "maf")
autoplot(loci_qc_Ivom, type = "missing")
#remove loci with greater than XX% missing data and less than X% MAF

Ivom_gt <- Ivom_gt %>% select_loci_if(loci_missingness(genotypes) < 0.1)
Ivom_gt <- Ivom_gt %>% select_loci_if(loci_maf(genotypes) > 0.05) # remove loci only found in a single copy in whole population
#Ivom_ig_gt <- Ivom_ig_gt %>% select_loci_if(loci_missingness(genotypes) < 0.1)
#Ivom_ig_gt <- Ivom_ig_gt %>% select_loci_if(loci_maf(genotypes) > 0.0015)
#update backing file
Ivom_gt <- gt_update_backingfile(Ivom_gt, backingfile = tempfile())
#Ivom_ig_gt <- gt_update_backingfile(Ivom_ig_gt, backingfile = tempfile())
# Impute missing data using mode
Ivom_gt_imp <- gt_impute_simple(Ivom_gt, method = "mode")
#Ivom_ig_gt_imp <- gt_impute_simple(Ivom_ig_gt, method = "mode")



###########################
## Population allele frequencies ##
###########################
Ivom_gt <- Ivom_gt %>%
  group_by(site)
pop_freqs <- loci_alt_freq(Ivom_gt)
Ivom_wild_gt <- Ivom_wild_gt %>%
  group_by(site)
pop_freqs <- loci_alt_freq(Ivom_wild_gt)

############
## PCA ##
############
# Ilex multispecies subset
redrep_list <- read_tsv("redrep.txt", col_names = FALSE)
redrep_gt <- Ilex_gt %>%
  filter(id %in% redrep_list$X1)
redrep_gt <- gt_update_backingfile(redrep_gt, backingfile = tempfile())
redrep_gt <- redrep_gt %>% select_loci_if(loci_maf(genotypes) > 0)
redrep_gt_imp <- gt_impute_simple(redrep_gt_imp, method = "mode")
redrep_gt_imp <- gt_update_backingfile(redrep_gt_imp, backingfile = tempfile())
gt_as_vcf(redrep_gt, "./Ilex_spp.vcf", overwrite = TRUE)

Ilex_pca <- gt_pca_partialSVD(redrep_gt_imp, k=7)
autoplot(Ilex_pca, type = "screeplot")
autoplot(Ilex_pca, type = "scores", k = c(1,2)) +
  aes(color = redrep_gt$spp) +
  labs(color = "spp")

pcs <- augment(x = Ilex_pca, data = redrep_gt)
eigenvalues <- tidy(Ilex_pca, "eigenvalues")

xlab <- paste("Axis 1 (", round(eigenvalues[1, 3], 1), " %)",
              sep = ""
)
ylab <- paste("Axis 2 (", round(eigenvalues[2, 3], 1), " %)",
              sep = ""
)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(8)

ggplot(data = pcs, aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = spp),
             shape = 21, size = 3, stroke = 0.1, show.legend = TRUE
  ) +
  #scale_fill_paletteer_d(palette="lisa::C_M_Coolidge") +
  scale_fill_manual(values = cols) +
  geom_text_repel(aes(label=id), size=1, max.overlaps=Inf, show.legend=FALSE) +
  labs(x = xlab, y = ylab)
#ggtitle("I. vomitoria PCA")

# Vomitoria only

Ivom_pca <- gt_pca_partialSVD(Ivom_gt_imp)
autoplot(Ivom_pca, type = "screeplot")
autoplot(Ivom_pca, type = "scores", k = c(1,2)) +
  aes(color = Ivom_gt$site) +
  labs(color = "site")

#intergenic PCA
Ivom_ig_gt <- Ivom_ig_gt %>% select_loci_if(loci_maf(genotypes) > 0)
Ivom_ig_pca <- gt_pca_partialSVD(Ivom_ig_gt)
autoplot(Ivom_ig_pca, type = "screeplot")
autoplot(Ivom_ig_pca, type = "scores", k = c(1,2)) +
  aes(color = Ivom_ig_gt$population) +
  labs(color = "superpopulation")
pcs <- augment(x = Ivom_ig_pca, data = Ivom_ig_gt)
write_tsv(pcs[, -2], "Ivom1-5_ig_pca.tsv")

#plot with ggplot
pcs <- augment(x = Ivom_pca, data = Ivom_gt)
eigenvalues <- tidy(Ivom_pca, "eigenvalues")
#eigenvalues_ig <- tidy(Ivom_ig_pca, "eigenvalues")
#pcs_ig <- augment(x = Ivom_ig_pca, data = Ivom_ig_gt)

xlab <- paste("Axis 1 (", round(eigenvalues[1, 3], 1), " %)",
              sep = ""
)
ylab <- paste("Axis 2 (", round(eigenvalues[2, 3], 1), " %)",
              sep = ""
)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(4)

ggplot(data = pcs, aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = population),
             shape = 21, size = 3, stroke = 0.1, show.legend = TRUE
  ) +
  #scale_fill_paletteer_d(palette="lisa::C_M_Coolidge") +
  scale_fill_manual(values = cols) +
  geom_text_repel(aes(label=id), size=1, max.overlaps=Inf, show.legend=FALSE) +
  labs(x = xlab, y = ylab)
  #ggtitle("I. vomitoria PCA")

# Calculate center for each population
centroid <- aggregate(cbind(.fittedPC1, .fittedPC2, .fittedPC3) ~ population,
                      data = pcs, FUN = mean
)

# Add these coordinates to our augmented pca object
pcs <- left_join(pcs, centroid, by = "population", suffix = c("", ".cen"))
ggplot(data = pcs, aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(xend = .fittedPC1.cen, yend = .fittedPC2.cen),
               show.legend = FALSE
  ) +
  geom_point(aes(fill = population),
             shape = 21, size = 3, show.legend = FALSE
  ) +
  geom_label(
    data = centroid,
    aes(label = population, fill = population),
    size = 4, show.legend = FALSE
  ) +
  scale_fill_manual(values = cols) +
  labs(x = xlab, y = ylab) +
  ggtitle("I. vomitoria PCA")

# add PC data to main dataset
Ivom_gt <- augment(Ivom_pca, data = Ivom_gt)
Ivom_gt %>% ggplot(aes(.fittedPC1, .fittedPC2, color = superpop)) +
  geom_point() +
  labs(x = "PC1", y = "PC2", color = "superpop")
#check loadings
autoplot(Ivom_pca, type = "loadings")

flagged_feral <- c("VA-FR-2","VA-FR-2X","FL-BE-1","VA-FR-1","VA-FR-1X","NC-CB-1","VA-FR-3","VA-FR-3X","NC-KW-2","NC-CB-2")

## non-cultivated IVV only ##
Ivom_wild_gt <- Ivom_gt %>%
  filter(population != "cultivated") %>%
  filter(id %!in% flagged_feral) %>%
  filter(site != "MC-AR")
Ivom_wild_gt <- Ivom_wild_gt %>% select_loci_if(loci_maf(genotypes) > 0)
Ivom_wild_gt <- gt_update_backingfile(Ivom_wild_gt, backingfile = tempfile())
Ivom_wild_gt_imp <- gt_impute_simple(Ivom_wild_gt, method = "mode")

Ivom_wild_gt$cluster <- wild_pops$cluster

Ivom_wild_pca <- gt_pca_partialSVD(Ivom_wild_gt_imp)
autoplot(Ivom_wild_pca, type = "screeplot")
autoplot(Ivom_wild_pca, type = "scores", k = c(1,3)) +
  aes(color = Ivom_wild_gt$site) +
  labs(color = "site")

pcs <- augment(x = Ivom_wild_pca, data = Ivom_wild_gt)
eigenvalues <- tidy(Ivom_wild_pca, "eigenvalues")

xlab <- paste("Axis 1 (", round(eigenvalues[1, 3], 1), " %)",
              sep = ""
)
ylab <- paste("Axis 2 (", round(eigenvalues[2, 3], 1), " %)",
              sep = ""
)

cols = gg_color_hue(3)

ggplot(data = pcs, aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = population),
             shape = 21, size = 3, stroke = 0.1, show.legend = TRUE
  ) +
  #scale_fill_paletteer_d(palette="lisa::C_M_Coolidge") +
  scale_fill_manual(values = cols) +
  #geom_text_repel(aes(label=id), size=1, max.overlaps=Inf, show.legend=FALSE) +
  labs(x = xlab, y = ylab)
#ggtitle("I. vomitoria PCA")

# subset out dune and inland inds
only_dune <- read_tsv("dune.txt", col_names = FALSE)
dune_gt <- Ivom_wild_gt %>%
  filter(id %in% only_dune$X1)
dune_gt <- dune_gt %>% select_loci_if(loci_maf(genotypes) > 0)

only_inland <- read_tsv("inland.txt", col_names = FALSE)
inland_gt <- Ivom_wild_gt %>%
  filter(id %in% only_inland$X1)
inland_gt <- inland_gt %>% select_loci_if(loci_maf(genotypes) > 0)
gt_as_vcf(dune_gt, "dune.vcf", overwrite = TRUE)
gt_as_vcf(inland_gt, "inland.vcf", overwrite = TRUE)

##############
## DAPC ##
##############

# visualize variance explained by PCs
tidy(Ivom_wild_pca, matrix = "eigenvalues") %>%
  ggplot(mapping = aes(x = PC, y = cumulative)) +
  geom_point()

Ivom_clusters <- gt_cluster_pca(Ivom_wild_pca, n_pca = 3, k_clusters = c(2,12))
autoplot(Ivom_clusters)
Ivom_clusters <- gt_cluster_pca_best_k(Ivom_clusters)
Ivom_clusters$best_k <- 3
Ivom_dapc <- gt_dapc(Ivom_clusters)
autoplot(Ivom_dapc, type = "scores")
autoplot(Ivom_dapc, type = "components", group = Ivom_gt$population)
autoplot(Ivom_dapc, "loadings")

# plot on map
Ivom_wild_gt <- Ivom_wild_gt %>% 
  ungroup() %>%
  mutate(dapc = Ivom_dapc$grp)
Ivom_wild_gt <- gt_add_sf(Ivom_wild_gt, c("lon", "lat"))
  
  
  ggplot() +
    geom_sf(data = map) +
    geom_sf(data = Ivom_wild_gt$geometry, aes(color = Ivom_wild_gt$dapc)) +
    coord_sf(
      xlim = c(-105, -75),
      ylim = c(20, 42)
    ) +
    labs(color = "DAPC cluster") +
    theme_minimal()

###############
## SNMF ##
###############
  
Ivom_snmf <- gt_snmf(
  x = Ivom_gt,
  k = 2:10,
  project = "new",
  n_runs = 3,
  entropy = TRUE,
  alpha = 50,
  percentage = 0.5,
  seed = c(239847,238973,568549)
)

autoplot(Ivom_snmf, type = "cv")
autoplot(Ivom_snmf, type = "barplot", k = 7, run = 1, annotate_group = TRUE,
         arrange_by_group = TRUE, arrange_by_indiv = TRUE,
         reorder_within_groups = TRUE)

q_mat <- get_q_matrix(Ivom_snmf, k = 7, run = 1) 
q_df <- q_mat %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = pops$ID,
         pop = pops$state,
         superpop = pops$superpop,
         order = pops$order)

q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with(".Q"), names_to = "prop", values_to = "q") 

q_df_prates <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))

q_df_prates %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = prop)) +
  #scale_fill_manual(values = q_palette, labels = c("AF", "Eam", "Wam")) +
  #scale_fill_viridis_d() +
  labs(fill = "cluster") +
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
  )

#plot the map
na_map <-
  ne_countries(continent = "north america", returnclass = "sf")

snmf_map <- ggplot(data = na_map) +
  geom_sf() +
  geom_point(
    data = anole_locs,
    aes(
      x = longitude,
      y = latitude,
      fill = population,
      size = 2
    ),
    shape = 21
  ) +
  scale_fill_viridis_d(direction = -1) +
  guides(size = FALSE,  #don't want a legend for the size
         fill = guide_legend(override.aes = list(size = 3))) +
  theme_bw(base_size = 16)
anole_map

#########################
## Heterozygosity ##
#########################
#Ivom_gt <- Ivom_wild_gt
Ivom_wild_gt <- Ivom_wild_gt %>% mutate(het_obs = indiv_het_obs(genotypes))
het <- Ivom_wild_gt %>% 
  group_by(site) %>%
  dplyr::summarize(mean = mean(het_obs))

write.csv(het, "./results_expanded_data/het_by_site_Ivom_wild.csv")

# add site coords
site_coords <- read.csv("coordinates_Ivom384.csv", header=TRUE)
het <- het %>%
  mutate(sitelats = site_coords[match(het$site,site_coords$Site),2],
       sitelons = site_coords[match(het$site,site_coords$Site),3]
  )
het <- st_as_sf(het, coords=c("sitelons", "sitelats"), crs = 4326)

ggplot() +
  geom_sf(data = map) +
  geom_sf(shape = 21, alpha = 0.8, data = het$geometry, aes(fill = het$mean)) +
  coord_sf(
    xlim = c(-105, -75),
    ylim = c(20, 42)
  ) +
  theme_minimal()

############################
## Mantel Test for IBD ##
############################

# generate allele sharing matrix
X <- attr(Ivom_wild_gt$genotypes, "fbm")
GSM <- snp_allele_sharing(X)

#visualize correlation
par(mar=c(4,4,0,0))
dens <- MASS::kde2d(Dgeodist, Dgen, n=300)
plot(Dgeodist, GSM, pch=20, cex=0.5,  
     xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(GSM ~ Dgeodist))
lines(loess.smooth(Dgeodist, GSM), col="red")

Dgen <- dist(GSM)
Dgeo <- dist(cbind(pcs$lat,pcs$lon), diag=T, upper=T)
Dgeodist <- geodist(cbind(pcs$lon,pcs$lat))
mtest <- vegan::mantel(Dgen,Dgeodist, method="spearman")

cor(Ivom_wild_gt$lon, Ivom_wild_gt$het_obs)

################
## Fst ##
################

# population Fst
single_sample <- c("FL-BW", "LA-JK","NC-KW","AR-TR","NC-CB")
Ivom_wild_gt <- Ivom_wild_gt %>% 
  group_by(site)
pop_Fst <- as_tibble(pop_fst(Ivom_wild_gt))
pop_Fst <- pop_Fst %>%
  filter(site %!in% single_sample)
write_tsv(pop_Fst, "pop_specific_Fst.tsv")
pop_Fst <- pop_Fst %>%
  mutate(site = het$site, geometry = het$geometry)

## pairwise Fst (Weir and Cockerham 1984) ##

pairwise_pop_Fst <- Ivom_wild_gt %>%
  group_by(site) %>%
  filter(site %!in% single_sample) %>%
  pairwise_pop_fst(type="tidy", method="WC84")

Ivom_wild_gt <- Ivom_wild_gt %>% group_by(cluster)
pairwise_cluster_Fst <- pairwise_pop_fst(Ivom_wild_gt, type="tidy", method="WC84")

write.csv(pairwise_cluster_Fst, "pairwise_pop_Fst_cluster_Ivomwild.csv")


pop_Fst <- pop_Fst %>%
  filter(value != "NaN")
pairwise_pop_Fst <- pairwise_pop_Fst %>%
  filter(value != "NaN")
pairwise_state_Fst <- pairwise_state_Fst %>%
  filter(value != "NaN")
single_ind_sites <- c("AR-TR","LA-JK","FL-BW","NC-CB")
pairwise_pop_Fst <- pairwise_pop_Fst %>%
  filter(site_1 %!in% single_ind_sites) %>%
  filter(site_2 %!in% single_ind_sites)
# plot just site Fst
ggplot() +
  geom_sf(data = map) +
  geom_sf(data = pop_Fst$geometry, aes(color = pop_Fst$value)) +
  scale_fill_viridis_c(option = "viridis", aesthetics = "color") +
  coord_sf(
    xlim = c(-105, -75),
    ylim = c(20, 42)
  ) +
  theme_minimal()

# pairwise Fst plotting
ggplot(data = pairwise_pop_Fst, aes(x=site_1, y=site_2, fill=value)) + 
  geom_tile() +
  #geom_text(data=pairwise_pop_Fst, aes(label=round(value, 2))) +
  theme_minimal() +
  scale_fill_gradient2(low = "white", high = "red", space = "Lab", name="Fst") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1))

# plot pairwise Fst with pop Fst
hist(pairwise_pop_Fst$value)
pairs_below_0.01 <- pairwise_pop_Fst %>% 
  filter(value < 0.01)

pairs_below_0.01 <- pairs_below_0.01 %>%
  mutate(site1lats = site_coords[match(pairs_below_0.01$site_1,site_coords$Site),2],
         site2lats = site_coords[match(pairs_below_0.01$site_2,site_coords$Site),2],
         site1lons = site_coords[match(pairs_below_0.01$site_1,site_coords$Site),3],
         site2lons = site_coords[match(pairs_below_0.01$site_2,site_coords$Site),3],
         pair_id = row.names(pairs_below_0.01)
         )
pairs_site1 <- pairs_below_0.01 %>%
  select(pair_id,site1lats,site1lons) %>%
  rename(lat = site1lats, lon = site1lons)
pairs_site2 <- pairs_below_0.01 %>%
  select(pair_id,site2lats,site2lons) %>%
  rename(lat = site2lats, lon = site2lons)

pairs <- bind_rows(pairs_site1, pairs_site2)

# plot sites colored by pop specific Fst and lines between pops with pairwise Fst < 0.01
ggplot() +
  geom_sf(data = map) +
  geom_sf(data = pop_Fst$geometry, aes(color = pop_Fst$value)) +
  scale_fill_viridis_c(option = "inferno", aesthetics = "color") +
  #geom_path(data=pairs, aes(x=lon, y=lat, group = pair_id), color="black", size=0.3, alpha = 0.7) +
  coord_sf(
    xlim = c(-105, -75),
    ylim = c(20, 42)
  ) +
  theme_minimal()

# calculate Fst in windows
Ivom_gt <- Ivom_gt %>% group_by(population)
window_pop_fst <- windows_pairwise_pop_fst(Ivom_gt, type = "tidy", 
                                           method = "WC84", window_size = 50000, 
                                           step_size = 10000, size_unit = "bp")


# get just atlantic-gulf comparisons
gulf_atl_fst <- window_pop_fst %>% 
  filter(stat_name == "fst_atlantic.gulf") %>%
  drop_na()

# get lengths of chromosomes and make cumulative base pair counts
fst_pos_cum <- gulf_atl_fst %>%
  group_by(chromosome) %>%
  summarise(max_bp = max(end)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(chromosome, bp_add)

# add cumulative counts back to original df
gulf_atl_fst <- gulf_atl_fst %>%
  inner_join(fst_pos_cum, by = "chromosome") %>%
  mutate(bp_cum = as.numeric(end) + as.numeric(bp_add)) 

# make list of chromosomes and their starting positions for x axis
fst_axis_set <- gulf_atl_fst %>%
  group_by(chromosome) %>%
  summarize(center = mean(bp_cum))

# get just atlantic-gulf comparisons
gulf_atl_fst <- window_pop_fst %>% 
  filter(stat_name == "fst_atlantic.gulf") %>%
  drop_na()

## do it all again for florida x gulf

# get just fl-gulf comparisons
fl_gulf_fst <- window_pop_fst %>% 
  filter(stat_name == "fst_florida.gulf") %>%
  drop_na()
# get lengths of chromosomes and make cumulative base pair counts
fst_pos_cum <- fl_gulf_fst %>%
  group_by(chromosome) %>%
  summarise(max_bp = max(end)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(chromosome, bp_add)

# add cumulative counts back to original df
fl_gulf_fst <- fl_gulf_fst %>%
  inner_join(fst_pos_cum, by = "chromosome") %>%
  mutate(bp_cum = as.numeric(end) + as.numeric(bp_add))

# make list of chromosomes and their starting positions for x axis
fst_axis_set <- fl_gulf_fst %>%
  group_by(chromosome) %>%
  summarize(center = mean(bp_cum))

# do it for all comparisons

# drop nas
fst_all <- window_pop_fst %>% 
  drop_na()
# get lengths of chromosomes and make cumulative base pair counts
fst_pos_cum <- fst_all %>%
  group_by(chromosome) %>%
  summarise(max_bp = max(end)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(chromosome, bp_add)

# add cumulative counts back to original df
fst_all <- fst_all %>%
  inner_join(fst_pos_cum, by = "chromosome") %>%
  mutate(bp_cum = as.numeric(end) + as.numeric(bp_add))

# make list of chromosomes and their starting positions for x axis
fst_axis_set <- fst_all %>%
  group_by(chromosome) %>%
  summarize(center = mean(bp_cum))

# the plot
genomewide_fst <- ggplot(fst_all, aes(
  x = bp_cum, y = value,
  color = as_factor(stat_name)
)) +
  facet_wrap(~stat_name) +
  geom_hline(
    yintercept = 0, color = "grey40",
    linetype = "dashed"
  ) +
  geom_point(alpha = 0.75) +
  scale_x_continuous(
    label = fst_axis_set$chromosome,
    breaks = fst_axis_set$center
  ) +
  labs(x = "Genomic Position", y = "FST") +
  #theme_minimal() + # for white background
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    #axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )
genomewide_fst

## plot by chromosome
fst_by_chrom <- gulf_atl_fst %>%
  filter(chromosome == "Chr02")

# the plot
chrom_fst <- ggplot(fst_by_chrom, aes(
  x = bp_cum, y = value,
  color = as_factor(stat_name)
)) +
  #  facet_wrap(~location) +
  geom_hline(
    yintercept = 0, color = "grey40",
    linetype = "dashed"
  ) +
  geom_line(alpha = 0.7) +
  scale_x_continuous(
    label = fst_axis_set$chromosome,
    breaks = fst_axis_set$center
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )
chrom_fst

###################################
## nucleotide diversity (pi) ##
###################################
species_pi <- Ivom_wild_gt %>%
  ungroup() %>%
  loci_pi()
mean(species_pi)

cluster_pi <- Ivom_wild_gt %>% 
  group_by(cluster) %>%
  loci_pi()

state_pi <- Ivom_gt %>%
  group_by(state) %>%
  loci_pi()

site_pi <- Ivom_wild_gt %>%
  group_by(site) %>%
  loci_pi()
# calculate averages
cluster_pi_avgs <- cluster_pi %>%
  group_by(group) %>%
  dplyr::summarize(avg_pi = mean(value))
state_pi_avgs <- state_pi %>%
  group_by(group) %>%
  summarize(avg_pi = mean(value))
site_pi_avgs <- site_pi %>%
  group_by(group) %>%
  dplyr::summarize(avg_pi = mean(value))

# map site pi
site_pi_avgs <- site_pi_avgs %>% 
  mutate(lat = site_coords[match(site_pi_avgs$group,site_coords$Site),2],
         lon = site_coords[match(site_pi_avgs$group,site_coords$Site),3]
         ) %>%
  rename(site = group)

site_pi_avgs <- st_as_sf(site_pi_avgs, coords = c(4,3), remove = FALSE, crs = 4326)
st_crs(map) == st_crs(site_pi_avgs)
site_pi_avgs <- site_pi_avgs %>%
  filter(avg_pi != "NaN")

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = site_pi_avgs$geometry, aes(color = site_pi_avgs$avg_pi)) +
  scale_fill_viridis_c(option = "inferno", aesthetics = "color") +
  coord_sf(
    xlim = c(-105, -75),
    ylim = c(20, 42)
  ) +
  theme_minimal()

#plot pi
pi_1 <- separate_wider_delim(superpop_pi, cols = loci, delim = "_", names = c("chr", "bp"))

data_cum <- pi_1 %>%
  group_by(chr) %>%
  summarise(max_bp = max(bp)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(chr, bp_add)

pi_1 <- pi_1 %>%
  inner_join(data_cum, by = "chr") %>%
  mutate(bp_cum = as.numeric(bp) + as.numeric(bp_add))

pi_axis_set <- pi_1 %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum))
#the plot
genomewide_pi <- ggplot(pi_1, aes(
  x = bp_cum, y = value,
  color = as_factor(group)
)) +
  facet_wrap(~group) +
  geom_hline(
    yintercept = 0, color = "grey40",
    linetype = "dashed"
  ) +
  geom_point(alpha = 0.75) +
  scale_x_continuous(
    label = pi_axis_set$chr,
    breaks = pi_axis_set$center
  ) +
  coord_cartesian(ylim = c(0, 1)) + 
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )
genomewide_pi

###################
## Tajima's D ##
###################

Ivom_noMAF_gt <- gen_tibble("Ilex1-5_noMAF.vcf.gz")
Ivom_noMAF_gt <- Ivom_noMAF_gt %>% left_join(pops, by = "id")
Ivom_noMAF_gt <- Ivom_noMAF_gt %>%
  filter(spp == "vomitoria") %>%
  filter(population != "cultivated") %>%
  filter(id %!in% flagged_feral) %>%
  filter(site != "MC-AR")
Ivom_noMAF_gt <- Ivom_noMAF_gt %>% select_loci_if(loci_maf(genotypes) > 0)
## estimate global Tajima's D values by superpop 
taj_d_by_site <- Ivom_noMAF_gt %>% 
  group_by(site) %>%
  pop_tajimas_d()

wild_sites <- Ivom_noMAF_gt$site %>%
  unique() %>%
  sort()

taj_d <- sapply(taj_d_by_site, `[[`, 1)

taj_d_by_site_2 <- tibble(wild_sites,taj_d)

taj_d_by_site_2 <- taj_d_by_site_2 %>% 
  mutate(lat = site_coords[match(taj_d_by_site_2$wild_sites,site_coords$Site),2],
         lon = site_coords[match(taj_d_by_site_2$wild_sites,site_coords$Site),3]
  )

taj_d_by_site_2 <- st_as_sf(taj_d_by_site_2, coords = c(4,3), remove = FALSE, crs = 4326)
st_crs(map) == st_crs(taj_d_by_site_2)
taj_d_by_site_2 <- taj_d_by_site_2 %>%
  filter(taj_d != "NaN")

taj_d_by_site_2 <- as_tibble(taj_d_by_site_2)

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = taj_d_by_site_2$geometry, aes(color = taj_d_by_site_2$taj_d)) +
  scale_fill_viridis_c(option = "inferno", aesthetics = "color") +
  coord_sf(
    xlim = c(-105, -75),
    ylim = c(20, 42)
  ) +
  theme_minimal()


# compute Tajima's D on a sliding window
sliding_TajD <- Ivom_gt %>% 
  group_by(population) %>%
  windows_pop_tajimas_d(window_size=20, step_size=10, size_unit="snp", min_loci=3)
# in bp windows
sliding_TajD <- Ivom_gt %>% 
  group_by(population) %>%
  windows_pop_tajimas_d(window_size=100000, step_size=50000, size_unit="bp", min_loci=3)


sliding_TajD2 <- sliding_TajD %>%
  pivot_longer(c("atlantic", "cultivated","florida","gulf"),
               names_to = "location")

sliding_TajD2 %>%
  ggplot(aes(x=start)) + 
  geom_line(aes(y=value, color=location), alpha=0.4) +
  facet_wrap(~chromosome)

#plotting genomewide data
pos_cum <- sliding_TajD2 %>%
  group_by(chromosome) %>%
  summarise(max_bp = max(end)) %>%
  mutate(snp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(chromosome, snp_add)

sliding_TajD3 <- sliding_TajD2 %>%
  inner_join(pos_cum, by = "chromosome") %>%
  mutate(snp_cum = end + snp_add)

axis_set <- sliding_TajD3 %>%
  group_by(chromosome) %>%
  summarize(center = mean(snp_cum))


#the plot
genomewide_TajD <- ggplot(sliding_TajD3, aes(
  x = snp_cum, y = value,
  color = as_factor(location)
)) +
  facet_wrap(~location) +
  geom_hline(
    yintercept = 0, color = "grey40",
    linetype = "dashed"
  ) +
  geom_line(alpha = 0.75) +
  scale_x_continuous(
    label = axis_set$chromosome,
    breaks = axis_set$center
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )
genomewide_TajD

#plot by chromosome
TajD_by_chrom <- sliding_TajD3 %>%
  filter(chromosome == "Chr01")

fst_by_chrom$value <- (scale(fst_by_chrom$value))
range(fst_by_chrom$value)
# scale data so they can be plotted together
TajD_by_chrom$value <- (scale(TajD_by_chrom$value))^2
range(TajD_by_chrom$value)


# the plot
chrom_TajD <- ggplot(data = TajD_by_chrom, aes(
  x = snp_cum, y = value
  #color = as_factor(location)
)) +
#  facet_wrap(~location) +
  geom_hline(
    yintercept = 0, color = "grey40",
    linetype = "dashed"
  ) +
  geom_line(alpha = 0.7) +
  scale_x_continuous(
    label = axis_set$chromosome,
    breaks = axis_set$center
  ) +
  #geom_line(data = fst_by_chrom, aes(x = bp_cum, y = value, color = "red")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )
chrom_TajD

########################
## Private alleles ##
########################

#get frequencies by X and turn them into binary data

allele_freqs <- Ivom_wild_gt %>%
  group_by(cluster) %>%
  loci_maf() %>%
  mutate(present = case_when(value > 0 ~ 1,
                             value == 0 ~ 0), 
        value = NULL
         )


#pivot data 
allele_freqs <- allele_freqs %>% pivot_wider(names_from = group, values_from = present)

#get alleles present only in one pop
private_alleles <- allele_freqs %>% 
  rowwise() %>%
  filter(sum(c_across(cols = 2:length(allele_freqs))) == 1)

private_by_site <- data.frame(colSums(private_alleles[,-1]))
write.csv(private_by_site, "private_by_site.csv")

private_by_pop <- data.frame(colSums(private_alleles[,-1]))
write.csv(private_by_pop, "private_by_pop.csv")

private_by_state <- data.frame(colSums(private_alleles[,-1]))
write.csv(private_by_state, "private_by_state.csv")

# plot on map
private_by_site <- private_by_site %>% 
  mutate(site = row.names(private_by_site))
private_by_site <- private_by_site %>% 
  mutate(lat = site_coords[match(private_by_site$site,site_coords$Site),2],
         lon = site_coords[match(private_by_site$site,site_coords$Site),3]
  ) %>%
  rename(private_allele_count = colSums.private_alleles....1..)

private_by_site <- st_as_sf(private_by_site, coords = c(4,3), remove = FALSE, crs = 4326)

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = private_by_site$geometry, aes(color = private_by_site$private_allele_count)) +
  scale_fill_viridis_c(option = "inferno", aesthetics = "color") +
  coord_sf(
    xlim = c(-105, -75),
    ylim = c(20, 42)
  ) +
  theme_minimal()

# plot by state
all_states <- map_data("state")
private_by_state <- private_by_state %>%
  mutate(region = tolower(state.name[match(rownames(private_by_state), state.abb)]))
  #rename(private_allele_count = colSums.private_alleles....1..)


ggplot(private_by_state, aes(map_id = region)) + 
  geom_map(aes(fill = private_allele_count), map = all_states) +
  scale_fill_viridis_c(option = "inferno", aesthetics = "color") +
  coord_sf(
    xlim = c(-105, -75),
    ylim = c(20, 42)
  ) +
  theme_minimal()

###########################
## LD decay plotting script 
###########################

rm(list = ls())
install.packages("tidyverse")
library(tidyverse)

# set path
my_bins <- "./Ivom_chr2.ld_decay_bins"

# read in data
ld_bins <- read_table(my_bins)

# plot LD decay
ggplot(ld_bins, aes(distance, avg_R2)) + geom_line() +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2)) +
  coord_cartesian(xlim = c(0,50000))

####################################################################
## RANGE EXPANSION ANALYSIS (Peter and Slatkin 2015, Evolution) ##
####################################################################

BiocManager::install("snpStats")
library(snpStats)
devtools::install_github("BenjaminPeter/rangeExpansion", ref="package")
library(rangeExpansion)

# import data
snp.file <- "Ivom_wild_rp_out_all.bed"
coord.file <- "rangeExp_coords_all.tsv"
region <- list(NULL, "REGION_1", "REGION_2", c("REGION_1", "REGION_2"))

count(rp_gt$site)

raw.data <- load.data(snp.file, n.snp=40000, coord.file, ploidy=2)
pop <- make.pop(raw.data, ploidy=2)
psi <- get.all.psi(pop)
psi_t <- t(psi)
psi_t[25 ,]
psi_t <- psi_t[-25,-25]
results <- run.single.region(region=c("REGION_1", "REGION_2"), pop=pop, psi=psi_t, xlen=20, ylen=40)
results_reg1 <- run.single.region(region="REGION_1", pop=pop, psi=psi_t, xlen=40, ylen=40)
results_reg2 <- run.single.region(region="REGION_2", pop=pop, psi=psi_t, xlen=10, ylen=20)

summary.origin.results(results)
plot.origin.results(results, add.samples = F)
summary.origin.results(results_reg1)
plot.origin.results(results_reg1)
summary.origin.results(results_reg2)
plot.origin.results(results_reg2)
colnames(psi) <- results[["coords"]][["pop"]]
row.names(psi) <- results[["coords"]][["pop"]]
pop_psi <- rowSums(psi)
pop_psi <- tibble(results[["coords"]][["pop"]],pop_psi)
colnames(pop_psi) <- c("site","psi")

pop_psi <- pop_psi %>%
  mutate(sitelats = site_coords[match(pop_psi$site,site_coords$Site),2],
         sitelons = site_coords[match(pop_psi$site,site_coords$Site),3]
  )
pop_psi <- st_as_sf(pop_psi, coords=c("sitelons", "sitelats"), crs = 4326)

pop_psi <- pop_psi %>%
  filter(site != "LA-JK")

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = pop_psi$geometry, aes(color = pop_psi$psi)) +
  scale_fill_viridis_c(option = "inferno", aesthetics = "color") +
  coord_sf(
    xlim = c(-105, -75),
    ylim = c(20, 42)
  ) +
  theme_minimal()

###################
## iHS-xpEHH ##
###################

install.packages("rehh")
library(rehh)
install.packages('R.utils')


# set chromosome
chrom <- "Chr01"
## read in data for each pop
# pop 1
atl_hh <- data2haplohh(hap_file = "dune_atl_phased.vcf.gz",
                         polarize_vcf = FALSE, chr.name = chrom,
                       )
# pop 2
gulf_hh <- data2haplohh(hap_file = "inland_atl_phased.vcf.gz",
                       polarize_vcf = FALSE, chr.name = chrom,
                       )

# filter on MAF
atl_hh_f <- base::subset(atl_hh, min_maf = 0.001)
gulf_hh_f <- subset(gulf_hh, min_maf = 0.001)

# perform scans
atl_scan <- scan_hh(atl_hh_f, polarized = FALSE)
gulf_scan <- scan_hh(gulf_hh_f, polarized = FALSE)

# compute EHH at specific location
#ehh <- calc_ehh(atl_hh_f, mrk=7591)
#plot(ehh)

# perform iHS on atlantic
atl_ihs <- ihh2ihs(atl_scan, freqbin = 1)
ggplot(atl_ihs$ihs, aes(POSITION, IHS)) + geom_point()
ggplot(atl_ihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point()

# perform xp-ehh
atl_gulf <- ies2xpehh(atl_scan, gulf_scan,
                       popname1 = "dune", popname2 = "inland"
                      )
ggplot(atl_gulf, aes(POSITION, XPEHH_dune_inland)) + geom_line()

xpehh_p <- quantile(atl_gulf$LOGPVALUE, 0.995, na.rm = TRUE)
ihs_atl_p <- quantile(atl_ihs$ihs[4], 0.995, na.rm = TRUE)
ihs_gulf_p <- quantile(gulf_ihs$ihs[4], 0.995, na.rm = TRUE)

cand_regions_ihs <- calc_candidate_regions(atl_ihs$ihs, threshold = ihs_p,
                                       window_size = 250000, overlap = 10000,
                                       min_n_extr_mrk = 3, pval = TRUE, ignore_sign = TRUE) %>%
                    select(CHR, START, END)

cand_regions_xp <- calc_candidate_regions(atl_gulf, threshold = xpehh_p,
                       window_size = 250000, overlap = 10000,
                       min_n_extr_mrk = 3, pval = TRUE, ignore_sign = TRUE) %>%
                    select(CHR, START, END)

cand_regions_overlap <- inner_join(cand_regions_ihs, cand_regions_xp)

distribplot(atl_gulf$XPEHH_atlantic_gulf, qqplot = TRUE)
manhattanplot(atl_gulf,
              main = "iHS (Atlantic x Gulf Population)")
manhattanplot(atl_gulf,
              main = paste("iHS ", chrom, "(Atlantic x Gulf Population)"),
              pval = TRUE)
# find the highest hits
hits <- atl_gulf %>% arrange(desc(LOGPVALUE)) %>% top_n(1)

# get SNP position
x <- hits$POSITION
x <- 47105349 # low xp-EHH

marker_id_a <- which(atl_hh@positions == x)
marker_id_g <- which(gulf_hh@positions == x)

# plot furcation patterns
atl_furcation <- calc_furcation(atl_hh, mrk = marker_id_a)
gulf_furcation <- calc_furcation(gulf_hh, mrk = marker_id_g)

plot(atl_furcation)
plot(gulf_furcation)

atl_haplen <- calc_haplen(atl_furcation)
gulf_haplen <- calc_haplen(gulf_furcation)

plot(atl_haplen)
plot(gulf_haplen)

# export xpEHH results
atl_gulf <- as_tibble(atl_gulf)
colnames(atl_gulf) <- tolower(colnames(atl_gulf))
file_name <- paste("./atl_gulf_xpEHH_", chrom, ".tsv", col="", sep="")
write_tsv(atl_gulf, file_name)

########################################
## Searching for candidate genes ##
########################################

# read in the gff
gff <-  read_table("Ilex_Hap1.filter.gff3", col_names = FALSE, comment = "#")
# subset and clear up the gff - add names
colnames(gff) <- c("chr", "source", "feature", "start", "end", "score",
                   "strand", "frame", "attribute")

# select genes only
gene_gff <- gff %>% 
  filter(feature == "gene")
# arrange the gff
gene_gff <- gene_gff %>% 
  arrange(start, end)
# make a gene mid point variable
gene_gff <- gene_gff %>% 
  mutate(mid = start + (end-start)/2)

# plot selection scan again
ggplot(atl_gulf, aes(position, logpvalue)) + geom_point()


# identify the highest peak of selection
hits <- atl_gulf %>% 
  arrange(desc(logpvalue)) %>% top_n(10)
# find highest hit
x <- hits$position[1]
# find hits closest to genes
gene_gff <- gene_gff %>% 
  mutate(hit_dist = abs(mid - x)) %>% 
  arrange(hit_dist)
# find hits within 250 Kb
gene_hits <- gene_gff %>% 
  mutate(hit_dist = abs(mid - x)) %>% 
  arrange(hit_dist) %>% 
  filter(hit_dist < 250000) %>%
  filter(chr == chrom)
# what are these genes?
gene_hits <- gene_hits %>% select(chr, start,end, attribute,hit_dist)
# separate out the attribute column
cand_names <- gene_hits %>% pull(attribute)

## by region
cand_regions <- cand_regions %>% 
  mutate(mid = START + (END-START)/2)
# find genes in region
gene_hits <- gene_gff %>% 
  filter(start > cand_regions$START - 1000, end < cand_regions$END + 1000) %>%
  filter(chr == chrom)
cand_names <- gene_hits %>% pull(attribute)


# get subsetted GFF
candidateGFF <- gff %>% 
  filter(attribute %in% cand_names)
gff_name <- paste("candidateGFF_", chrom, ".gff3", col="", sep="")
write_tsv(candidateGFF, gff_name, col_names = FALSE)

############################################################
## loop for running iPS and XP-EHH by chromosome ##
############################################################

chrlist <- c()
n <- 20
ihsatllist <- vector("list", length = n)
ihsgulflist <- vector("list", length = n)
xplist <- vector("list", length = n)
rsblist <- vector("list", length = n)
# make list of chromosome names
for (i in 1:n) {
  chrom <- ifelse(i < 10, paste("Chr0",i, sep=""), paste("Chr",i, sep=""))
  chrlist[i] <- chrom
}

# iPS, XP-EHH, Rsb loop
for (i in 1:n){
  # pop 1
  atl_hh <- data2haplohh(hap_file = "dune_gulf_phased.vcf.gz",
                         polarize_vcf = FALSE, chr.name = chrlist[i],
  )
  # pop 2
  gulf_hh <- data2haplohh(hap_file = "inland_gulf_phased.vcf.gz",
                          polarize_vcf = FALSE, chr.name = chrlist[i],
  )
  
  # filter on MAF
  atl_hh_f <- subset(atl_hh, min_maf = 0.001)
  gulf_hh_f <- subset(gulf_hh, min_maf = 0.001)
  
  # perform scans
  atl_scan <- scan_hh(atl_hh_f, polarized = FALSE)
  gulf_scan <- scan_hh(gulf_hh_f, polarized = FALSE)
  
  # perform iHS 
  atl_ihs <- ihh2ihs(atl_scan, freqbin = 1)
  gulf_ihs <- ihh2ihs(gulf_scan, freqbin = 1)
  
  # perform xp-ehh
  atl_gulf <- ies2xpehh(atl_scan, gulf_scan,
                        popname1 = "Atlantic", popname2 = "Gulf")
  # perform Rsb
  rsb_atl_gulf <- ines2rsb(scan_pop1 = atl_scan,
                           scan_pop2 = gulf_scan,
                           popname1 = "Atlantic",
                           popname2 = "Gulf")
  
  ihsatllist[[i]] <- atl_ihs$ihs
  ihsgulflist[[i]] <- gulf_ihs$ihs
  xplist[[i]] <- atl_gulf
  rsblist[[i]] <- rsb_atl_gulf
}

ihs_atl_wg <- do.call(rbind, ihsatllist)
ihs_gulf_wg <- do.call(rbind, ihsgulflist)
xpehh_wg <- do.call(rbind, xplist)
rsb_wg <- do.call(rbind, rsblist)

xpehh_p <- quantile(xpehh_wg$LOGPVALUE, 0.9995, na.rm = TRUE)
ihs_atl_p <- quantile(ihs_atl_wg[4], 0.9995, na.rm = TRUE)
ihs_gulf_p <- quantile(ihs_gulf_wg[4], 0.9995, na.rm = TRUE)
rsb_p <- quantile(rsb_wg[4], 0.9995, na.rm = TRUE)

colnames(xpehh_wg) <- c("CHR", "POSITION", "XPEHH_dune_inland", "LOGPVALUE")
manhattanplot(xpehh_wg,
              main = "XP-EHH (Gulf Dune x Inland Populations)",
              pval = FALSE,
              threshold = c(100),
              cr = cand_regions_xp_gulf)
manhattanplot(ihs_atl_wg,
              main = "iHS (Dune Populations)",
              pval = TRUE,
              threshold = ihs_atl_p)
manhattanplot(ihs_gulf_wg,
              main = "iHS (Inland Populations)",
              pval = TRUE,
              threshold = ihs_gulf_p)
manhattanplot(rsb_wg,
              main = "Rsb (Atlantic x Gulf Population)",
              pval = TRUE,
              threshold = rsb_p)

cand_regions_ihs_atl <- calc_candidate_regions(ihs_atl_wg, threshold = ihs_atl_p,
                                           window_size = 250000, overlap = 10000,
                                           min_n_extr_mrk = 3, pval = TRUE, ignore_sign = TRUE) %>%
                                          dplyr::select(CHR, START, END)

cand_regions_ihs_gulf <- calc_candidate_regions(ihs_gulf_wg, threshold = ihs_gulf_p,
                                               window_size = 250000, overlap = 10000,
                                               min_n_extr_mrk = 3, pval = TRUE, ignore_sign = TRUE) %>%
                                              dplyr::select(CHR, START, END)

cand_regions_xp <- calc_candidate_regions(xpehh_wg, threshold = 7, negativeThreshold = -7,
                                          window_size = 250000, overlap = 10000,
                                          min_n_extr_mrk = 3, pval = FALSE, ignore_sign = FALSE)
                                          dplyr::select(CHR, START, END) 

cand_regions_rsb <-  calc_candidate_regions(rsb_wg, threshold = rsb_p,
                                            window_size = 250000, overlap = 10000,
                                            min_n_extr_mrk = 3, pval = TRUE, ignore_sign = TRUE) %>%
                                            dplyr::select(CHR, START, END)                         

# cand_regions_overlap_atl <- inner_join(cand_regions_ihs_atl, cand_regions_xp) %>%
#                             mutate(pop = "Atlantic")
# cand_regions_overlap_gulf <- inner_join(cand_regions_ihs_gulf, cand_regions_xp) %>%
#                               mutate(pop = "Gulf")
# cand_regions_overlap_rsb <- inner_join(cand_regions_rsb, cand_regions_xp) %>%
#   mutate(pop = "Rsb")
# 
# 
# cand_regions_overlap <- bind_rows(cand_regions_overlap_atl, cand_regions_overlap_gulf) %>%
#   distinct() %>%
#   rename(chromosome = CHR, start = START, end = END) %>%
#   arrange(chromosome, start)
# 
# gulf_atl_fst <- gulf_atl_fst %>%
#   arrange(chromosome, start)

# ihs, rsb, xp-ehh

cand_regions_ihs <- bind_rows(cand_regions_ihs_atl, cand_regions_ihs_gulf)
# 1. Find rows common to df1 and df2, but not df3
in_ihs_rsb <- inner_join(cand_regions_ihs, cand_regions_rsb) %>% anti_join(cand_regions_xp)

# 2. Find rows common to df1 and df3, but not df2
in_ihs_xp <- inner_join(cand_regions_ihs, cand_regions_xp) %>% anti_join(cand_regions_rsb)

# 3. Find rows common to df2 and df3, but not df1
in_rsb_xp <- inner_join(cand_regions_rsb, cand_regions_xp) %>% anti_join(cand_regions_ihs)

# Combine the results
cand_regions_overlap_2 <- bind_rows(in_ihs_rsb, in_ihs_xp, in_rsb_xp)

# 4. Find rows common to all three
in_ihs_rsb_xp <- inner_join(cand_regions_ihs, cand_regions_rsb) %>% inner_join(cand_regions_xp)

cand_regions_overlap_23 <- bind_rows(cand_regions_overlap_2, in_ihs_rsb_xp)

# find genes in region
genes_list <- vector("list", length = nrow(cand_regions_xp_atl))
for (i in 1:nrow(cand_regions_xp_atl)) {
  hits <- ann_gene_gff %>% 
    filter(chr == cand_regions_xp_atl[i,1], start > cand_regions_xp_atl[i,2]-1000, end < cand_regions_xp_atl[i,3]+1000)
  genes_list[[i]] <- hits
}

genes_list <- vector("list", length = nrow(cand_regions_xp_gulf))
for (i in 1:nrow(cand_regions_xp_gulf)) {
  hits <- ann_gene_gff %>% 
    filter(chr == cand_regions_xp_gulf[i,1], start > cand_regions_xp_gulf[i,2]-1000, end < cand_regions_xp_gulf[i,3]+1000)
  genes_list[[i]] <- hits
}
# for putting all regions together
gene_hits <- bind_rows(genes_list)

# store hits and regions
gene_hits_atl <- gene_hits
gene_hits_gulf <- gene_hits
cand_regions_xp_atl <- cand_regions_xp
cand_regions_xp_gulf <- cand_regions_xp

gene_hits_atl$ann <- substr(gene_hits_atl$ann, 6, 1000)
gene_hits_gulf$ann <- substr(gene_hits_gulf$ann, 6, 1000)

gene_hits_atl_nona <- gene_hits_atl %>%
  filter(ann != "NA")

gene_hits_overlap_atl_gulf <- inner_join(gene_hits_atl, gene_hits_gulf)
write_tsv(gene_hits_overlap_atl_gulf, "gene_hits_overlap_atl_gulf.tsv", col_names = FALSE)
write_tsv(gene_hits_atl_nona, "gene_hits_EHH_atl.tsv", col_names = FALSE)
write_tsv(gene_hits_gulf, "gene_hits_EHH_gulf.tsv", col_names = FALSE)

# make background gene set #
random_gene_gff <- ann_gene_gff[sample(nrow(ann_gene_gff), 5000), ]
random_gene_gff$ann <- substr(random_gene_gff$ann, 6, 1000)
write_tsv(random_gene_gff, "background_gene_set5000.tsv", col_names = FALSE)

# make regions file for querying using bcftools
cand_regions_xp_atlgulf <- bind_rows(cand_regions_xp_atl,cand_regions_xp_gulf)
cand_regions_xp$region <- NULL
write_tsv(cand_regions_xp_atlgulf, "cand_regions_xp_atlgulf.tsv", col_names = TRUE)
# find snps overlapping between EHH and RDA
EHH_snps <- read_tsv("cand_region_snps_xp_atlgulf.txt")
EHH_snps <- EHH_snps %>%
  select(CHROM,POS,ID)
RDA_snps <- read.csv("./results_expanded_data/Chapter2/rda_cand_snps.csv")
EHH_snps <- EHH_snps %>% mutate(snp = paste(CHROM,"_",POS, sep=""))
RDA_snps$snp <- substr(RDA_snps$snp,1,nchar(RDA_snps$snp)-2)
RDA_coast_snps <- RDA_snps %>% filter(predictor == "COAST")
EHH_RDA_overlap <- RDA_snps %>%
  filter(snp %in% EHH_snps$snp)



# another way
EHH_snps <- read_tsv("EHH_cand_regions_snps.tsv")

colnames(EHH_snps)[3] <- "snp"
cand$snp <- substr(cand$snp,1,nchar(cand$snp)-2)

RDA_EHH_overlap <- inner_join(EHH_snps, cand, by = "snp")

# extract candidate regions around SNPs
EHH_RDA_overlap <- separate_wider_delim(EHH_RDA_overlap, 
                     cols = snp, delim = "_", names = c("CHROM", "POS"))
EHH_RDA_overlap$POS <- as.numeric(EHH_RDA_overlap$POS)
overlap_regions <- EHH_RDA_overlap %>%
  mutate(START = POS - 12500,
         END = POS + 12500)

# find genes in region
genes_list <- vector("list", length = nrow(overlap_regions))
for (i in 1:nrow(overlap_regions)) {
  chrom <- as.character(overlap_regions[i,3])
  lower_bound <- as.numeric(overlap_regions[i,15]-1000)
  upper_bound <- as.numeric(overlap_regions[i,16]+1000)
  hits <- ann_gene_gff %>% 
    filter(chr == chrom, start > lower_bound, end < upper_bound)
  genes_list[[i]] <- hits
}  

# combine lists and remove duplicates
RDA_EHH_gene_hits <- bind_rows(genes_list)
length(RDA_EHH_gene_hits$attribute[duplicated(RDA_EHH_gene_hits$attribute)])
RDA_EHH_gene_hits <- RDA_EHH_gene_hits[!duplicated(RDA_EHH_gene_hits$attribute),]

RDA_EHH_gene_hits$ann <- substr(RDA_EHH_gene_hits$ann, 6, 1000)
write_tsv(RDA_EHH_gene_hits, "RDA_EHH_gene_hits.gff3", col_names = FALSE)  


#############
## test whether RDA SNPs are more elevated in EHH stats than random SNPs ##
#############
xpehh_wg <- xpehh_wg %>% mutate(snp = paste(CHR,"_",POSITION, sep=""))

RDA_in_xp <- xpehh_wg %>% filter(snp %in% RDA_snps$snp)
mean(RDA_in_xp$LOGPVALUE, na.rm = TRUE)
RDA_coast_in_xp <- xpehh_wg %>% filter(snp %in% RDA_coast_snps$snp)
mean(RDA_coast_in_xp$LOGPVALUE, na.rm = TRUE)

random_xp <- xpehh_wg[sample(nrow(xpehh_wg), 10000), ]
mean(random_xp$LOGPVALUE, na.rm = TRUE)

t.test(random_xp$XPEHH_Atlantic_Gulf, RDA_coast_in_xp$XPEHH_Atlantic_Gulf)

xpehh_p <- quantile(RDA_in_xp$LOGPVALUE, 0.9995, na.rm = TRUE)

manhattanplot(RDA_coast_in_xp,
              main = "XP-EHH (Atlantic x Gulf Population)",
              #pval = TRUE,
              threshold = xpehh_p)
manhattanplot(random_xp,
              main = "XP-EHH (Atlantic x Gulf Population)",
              #pval = TRUE,
              threshold = xpehh_p)

# extract candidate regions around SNPs
block <- 12500
coast_cand_snps <- RDA_coast_in_xp %>%
  filter(LOGPVALUE > 2)
coast_cand_regions <- coast_cand_snps %>%
  mutate(START = POSITION - 25000,
         END = POSITION + 25000)

RDA_coast_snps <- separate_wider_delim(RDA_coast_snps, 
                             cols = snp, delim = "_", names = c("chr", "pos"))

coast_RDA_regions <- RDA_coast_snps %>%
  mutate(START = as.numeric(pos) - block,
         END = as.numeric(pos) + block)

# find genes in region
genes_list <- vector("list", length = nrow(cand_regions_overlap_23))
for (i in 1:nrow(coast_cand_regions)) {
  hits <- gene_gff %>% 
    filter(chr == coast_cand_regions[i,1], start > coast_cand_regions[i,6]-1000, end < coast_cand_regions[i,7]+1000)
  genes_list[[i]] <- hits
}  

genes_list <- vector("list", length = nrow(coast_RDA_regions))
for (i in 1:nrow(coast_RDA_regions)) {
  hits <- gene_gff %>% 
    filter(chr == coast_RDA_regions[i,3], start > coast_RDA_regions[i,14], end < coast_RDA_regions[i,15])
  genes_list[[i]] <- hits
}  

# combine lists of genes and remove duplicates
coast_gene_hits <- bind_rows(genes_list)
length(coast_gene_hits$attribute[duplicated(coast_gene_hits$attribute)])
coast_gene_hits <- coast_gene_hits[!duplicated(coast_gene_hits$attribute),]
coast_RDA_gene_hits <- bind_rows(genes_list)

write_tsv(coast_gene_hits, "coast_gene_hits.gff3", col_names = FALSE)  

rsb_wg <- rsb_wg %>% mutate(snp = paste(CHR,"_",POSITION, sep=""))
RDA_in_rsb <- rsb_wg %>% filter(snp %in% RDA_snps$snp)
RDA_coast_in_rsb <- rsb_wg %>% filter(snp %in% RDA_coast_snps$snp)
RDA_in_either_rsb_xp <- bind_rows(RDA_coast_in_rsb,RDA_coast_in_xp)


gene_hits_chr2 <- genes_list[[6]]
gene_hits_chr3 <- genes_list[[1]]
gene_hits_chr4 <- genes_list[[2]]
gene_hits_chr6 <- genes_list[[3]]
gene_hits_chr12 <- genes_list[[5]]  
gene_hits_chr13 <- genes_list[[4]]




# for dealing with orthologs from OrthoFinder
###

gene_hits_chr3 <- separate_wider_delim(gene_hits_chr3, cols = attribute, delim = ";", names = c("ID", "locus"))

gene_hits_chr3$ID <- substr(gene_hits_chr3$ID, 4, 16)
cand_orthologs <- orthologs %>%
  filter(grepl(Ilex_Hap1_peptides.filter, gene_hits_chr3$ID))

cand_orthologs <- subset(
  orthologs,
  grepl(
    paste0(gene_hits_chr3$ID, collapse = "|"),
    orthologs$Ilex_Hap1_peptides.filter,
    ignore.case = TRUE
  )
)

cand_orthologs <- cand_orthologs %>%
  separate_wider_delim(cols = Hannuus_494_r1.2.protein, names = c("first","secondary"), 
                       delim = ",", too_many = "drop", too_few = "align_start", cols_remove = FALSE)

write.csv(cand_orthologs, "cand_orthologs_chr3.csv")

###


cand_names <- gene_hits_chr2 %>% 
  pull(attribute) 

# get subsetted GFF
candidateGFF <- gff %>% 
  filter(attribute %in% cand_names) %>%
  distinct()
gff_name <- "gulfxatl_chr2.gff3"
write_tsv(candidateGFF, gff_name, col_names = FALSE)

## Venn diagrams
install.packages("ggVennDiagram")
library(ggVennDiagram)
library(purrr)
df_list <- list(cand_regions_ihs,cand_regions_rsb,cand_regions_xp)
purrr::reduce(df_list, inner_join)
overlap <- calculate.overlap(
  x = list(
    "iHS" = cand_regions_ihs, 
    "Rsb" = cand_regions_rsb, 
    "XP-EHH" = cand_regions_xp)
)
names(overlap) <- c("a123", "a12", "a13", "a23", "a1", "a2", "a3")

ggVennDiagram(x=list(cand_regions_ihs$START, cand_regions_rsb$START, cand_regions_xp$START),
              category.names = c("iHS","Rsb","XP-EHH"))

## get single copy orthologs from Orthofinder results
orthologs_Ipa <- read_tsv("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1-5/CAUOFW02_translation__v__Ilex_Hap1_peptides.filter.tsv"
                     )
orthologs <- separate_wider_delim(orthologs, cols = CAUOFW02_translation, delim = ",", names_sep = "_", too_few = "align_start")
orthologs <- separate_wider_delim(orthologs, cols = Ilex_Hap1_peptides.filter, delim = ",", names_sep = "_", too_few = "align_start")
sc_orthologs <- orthologs %>%
  filter(is.na(CAUOFW02_translation_2)) %>%
  filter(is.na(Ilex_Hap1_peptides.filter_2))
sc_orthologs_list <- as_tibble(sc_orthologs$Ilex_Hap1_peptides.filter_1)
write_tsv(sc_orthologs_list, "sc_orthologs_Ivo.txt", col_names = FALSE)

sc_orthologs_list$value <- substr(sc_orthologs_list$value, 2, 13)

# subset gff using list of gene ids
ortho_pattern <- paste(sc_orthologs_list$value, collapse = "|")

sc_orthologs_gff <- gff %>% 
  filter(str_detect(attribute, ortho_pattern)) 

write_tsv(sc_orthologs_gff, "Ivo_Ipa_singlecopy.gff3", col_names = FALSE)

## Gene ontology
# get reference set
ref_set <- gene_gff[sample(nrow(gene_gff), 10000), ]
write_tsv(ref_set, "ref_set_gff.gff3", col_names = FALSE)

###################
## M-K Test ##
###################
install.packages("jsonlite")
install.packages("janitor")
library(jsonlite)
library(janitor)
mk_imputed <- fromJSON(txt ="./mkado_imputed_43843288.out",
         simplifyDataFrame = TRUE,
         flatten = FALSE) %>%
  as_tibble() %>%
  t()
gene <- row.names(mk_imputed)
mk_imputed <- as_tibble(mk_imputed)
mk_imputed <- tibble(gene, mk_imputed)
colnames(mk_imputed) <- c("gene","alpha","p_val","pn_neutral","pwd","dn","ds","pn_total","ps_total","cutoff")
mk_imputed$p_val <- as.numeric(mk_imputed$p_val)

# from I. dis (unfiltered dataset)
mk_imputed_Idis <- fromJSON(txt ="./mkado_imputed_Idis_43911533.out",
                       simplifyDataFrame = TRUE,
                       flatten = FALSE) %>%
  as_tibble() %>%
  t()
gene <- row.names(mk_imputed_Idis)
mk_imputed_Idis <- as_tibble(mk_imputed_Idis)
mk_imputed_Idis <- tibble(gene, mk_imputed_Idis)
colnames(mk_imputed_Idis) <- c("gene","alpha","p_val","pn_neutral","pwd","dn","ds","pn_total","ps_total","cutoff")
mk_imputed_Idis$p_val <- as.numeric(mk_imputed_Idis$p_val)

# vanilla mk test from degenotate
mk <- read_tsv("./degenotate/mk.tsv")
mk <- mk %>%
  arrange(pval) %>%
  mutate(top = FALSE)
mk_cand <- head(mk, 20)
mk <- mk %>% 
  mutate(top = case_when(
    transcript %in% mk_cand$transcript ~ TRUE,
    transcript %!in% mk_cand$transcript ~ FALSE)
        )

mk_high_a <- mk %>%
  filter(CI_low > 0.7)

mk_no_na <- mk %>%
  drop_na() 
mk_top_20 <- mk %>%
  arrange(pval) %>%
  head(20)

mean(mk_cand$pN)
mean(mk_cand$pS)
mean(mk_cand$dN)
mean(mk_cand$dS)

fixed <- c(0,4)
poly <- c(7,20)
cont <- tibble(fixed, poly)
fisher.test(cont)

## conduct Benjamini-Hochberg correction ##
hist(mk_imputed$p_val)
hist(-log(mk_no_na$pval))
hist(mk_no_na$dos)
hist(mk_no_na$odds_ni)

# conduct B-H correction and extract significant values
mk <- mk %>%
  mutate(p_adj = p.adjust(pval, method = "fdr"))


mk_imputed_cands <- mk_imputed %>%
  filter(p_adj < 0.10)

mk_imputed <- mk_imputed %>%
  mutate(p_adj = p.adjust(p_val, method = "fdr"))
mk_imputed_cands <- mk_imputed %>%
  filter(p_adj < 0.10)

mk_imputed_Idis <- mk_imputed_Idis %>%
  mutate(p_adj = p.adjust(p_val, method = "fdr"))
mk_imputed_Idis_cands <- mk_imputed_Idis %>%
  filter(p_adj < 0.10)

# vanilla test from mkado
mkado <- read_tsv("./mkado_results_Ipa_normal.tsv")

volcano_mk <- mk %>% #create plot
  ggplot(aes(x = dos, y = -log10(pval), label = transcript, col = top)) +
  #geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  #geom_hline(yintercept = -log10(0.001), col = "gray", linetype = 'dashed') +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  guides(col="none")
  #scale_color_manual(values = c("forestgreen", "tan3", "grey"), # to set the colours of our variable
  #                   labels = c("Down in roasted", "Up in roasted", "Not significant")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  #coord_cartesian(ylim = c(0, 20), xlim = c(-1, 1)) + # since some genes can have minuslog10padj of inf, we set these limits
  #labs(color = '', #legend_title
  #     x = expression("alpha"), y = expression("-log"[10]*"p-value")) +
  #scale_x_continuous(breaks = seq(-3, 7, 1)) + # to customise the breaks in the x axis
  #ggtitle('Chemicals in negative mode significantly affected by the roasting process') + # Plot title
  #geom_text_repel(size = 6, max.overlaps = 10, show.legend = FALSE)  # To show all labels 
  #vol_theme +
  #theme(legend.position="none")
volcano_mk 

mk_imputed_cands$gene <- substr(mk_imputed_cands$gene, 2, 13)
mk_imputed_GFF <- gene_gff %>% 
  filter(ID %in% mk_imputed_cands$gene) %>%
  distinct()
gff_name <- "mk_imputed_cand.gff3"
write_tsv(mk_imputed_GFF, gff_name, col_names = FALSE)

gene_gff <- separate_wider_delim(gene_gff, cols = attribute, delim = ";", names = c("ID", "locus"))
mk_top_20$transcript <- substr(mk_top_20$transcript, 2, 13)
ann_gene_gff$ID <- substr(ann_gene_gff$ID, 2, 13)

mkGFF <- ann_gene_gff %>% 
  filter(ID %in% mk_top_20$transcript) %>%
  distinct()
gff_name <- "mk_top20.gff3"
mkGFF <- mkGFF %>%
  mutate(ann = ann$sseqid)
write_tsv(mkGFF, gff_name, col_names = FALSE)

write_tsv(tibble(gene_gff$ID), "gene_IDs.txt", col_names = FALSE)

## Add functional annotations from DIAMOND BLAST ##

blast_results <- read_tsv("gene_annotations_swissprot.txt", col_names = FALSE)
colnames(blast_results) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast_results_1 <- blast_results %>% 
  distinct(qseqid, .keep_all = TRUE)
blast_results_1 %>%
  select(sseqid) %>%
  write_tsv("ann_gis.txt")
gene_gff <- gff %>% 
  filter(feature == "gene")
gene_gff <- gene_gff %>% 
  arrange(start, end)
gene_gff <- separate_wider_delim(gene_gff, cols = attribute, delim = ";", names = c("ID", "locus"))
gene_gff$ID <- substr(gene_gff$ID, 4, 16)
blast_results_1$qseqid <- substr(blast_results_1$qseqid, 1, 13)
ann_gene_gff <- gene_gff %>% 
  mutate(ann = blast_results_1[match(gene_gff$ID,blast_results_1$qseqid),2])
ann_gene_gff <- ann_gene_gff %>%
  mutate(ann = paste("Name=",ann$sseqid, sep = ""))
ann_gene_gff <- ann_gene_gff %>%
  mutate(attribute = paste(ID,locus,ann, sep=";"))

