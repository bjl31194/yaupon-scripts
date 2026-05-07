install.packages("raster")
install.packages("sp")
install.packages("geodata")
install.packages("psych")
install.packages("ggOceanMaps")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("mapmixture")
install.packages("gridExtra")
install.packages("devtools")
install.packages("rnaturalearthhires")
install.packages("RColorBrewer")
install.packages("maps")
install.packages("mapplots")
BiocManager::install("qvalue")
install.packages("geometry")
install_github("landscape-genomics/rdadapt")
library(rdadapt)
library(RColorBrewer)
library(mapmixture)
library(geodata)
library(raster)
library(sp)
library(terra)
library(psych)
library(vegan)
library(tidyverse)
library(tibble)
library(dplyr)
library(plyr)
library(ggOceanMaps)
library(rnaturalearth)
library(rnaturalearthdata)
library(gridExtra)
library(ggplot2)
library(devtools)
library(rnaturalearthhires)
library(maps)
library(mapplots)

Americas <- ne_countries(scale = 50, continent = c("North America","South America"))
medium_scale_map <- ggplot() +
  geom_sf(data = Americas) +
  coord_sf(xlim = c(-140, -25), ylim = c(-40, 40))
medium_scale_map

#########################################
## Plot structure results onto map ##
#########################################

# Create admixture file format 2
admixture <- read.csv("admix_Ivom384_2.csv")


admixture2 <- admixture %>% 
  mutate(
    max = apply(admixture[3:6],1,max),
    Assignment = colnames(admixture)[apply(admixture,1,which.max)],
    Admix = case_when(
      max <= 0.6 ~ TRUE,
      TRUE ~ FALSE,
    )
  )

admixture2 %>%
  mutate(Assignment = replace(Assignment, Admix == TRUE, paste(Assignment, 'weak')))

require(tidyr)
admixture_by_ind <- admixture2 %>%
  mutate(Assignment = replace(Assignment, Admix == TRUE, paste(Assignment, 'weak', sep = ''))) %>%
  group_by(Site, Assignment) %>%
  summarize(prop = n()) %>%
  mutate(prop = prop/sum(prop)) %>%
  pivot_wider(names_from = Assignment,  values_from = prop) %>%
  replace(is.na(.),0) %>%
  mutate(Ind = Site, .after = Site) %>%
  relocate(c(Cluster2,Cluster3,Cluster4,Cluster1weak,Cluster2weak,Cluster3weak,Cluster4weak), .after = Cluster1)

# Read in admixture file format 1
file1 <- system.file("extdata", "admixture1.csv", package = "mapmixture")
admixture1 <- read.csv("admix_Ivom384.csv")
admix_example <- read.csv(file1)
# Read in coordinates file
file2 <- system.file("extdata", "coordinates.csv", package = "mapmixture")
coordinates <- read.csv("coordinates_Ivom384.csv")
coord_example <- read.csv(file2)

coordinates <- coordinates[-4]

#filter out sites
sites_to_remove <- c("MCIL","FL-7")
coordinates <- coordinates[!(coordinates$Site %in% sites_to_remove), ]
admixture1 <- admixture1[!(admixture1$Site %in% sites_to_remove), ]

# Run mapmixture
map1 <- mapmixture(admixture_offset, coordinates1, crs = 3035)
map1

# Run mapmixture (chad version, by mean membership per pop)
target <- c("United States of America","mexico")
map2 <- mapmixture(
  admixture_df = admixture1,
  coords_df = coordinates,
  cluster_cols = brewer.pal(4, "Set1"),
  cluster_names = c("Cluster1","Cluster2","Cluster3","Cluster4"),
  crs = 4326,
  boundary = c(xmin=-98, xmax=-73, ymin=26, ymax=37),
  #boundary = c(xmin=-90, xmax=-75, ymin=25, ymax=35),
  pie_size = 0.5,
  pie_border = 0.5,
  pie_border_col = "black",
  pie_opacity = 0.9,
  land_colour = "#d9d9d9",
  sea_colour = "#deebf7",
  basemap = rnaturalearth::ne_states(country=c("United States of America","mexico")),
  expand = TRUE,
  arrow = TRUE,
  arrow_size = 1.5,
  arrow_position = "bl",
  scalebar = TRUE,
  scalebar_size = 1.5,
  scalebar_position = "tl",
  plot_title = "Admixture Map",
  plot_title_size = 12,
  axis_title_size = 10,
  axis_text_size = 8
)
map2 

# offset pies that overlap
coordinates1 <- coordinates[16 ,]
coordinates2 <- coordinates1 %>% 
  mutate(Lat1 = 28.5, Lon1 = -83.5)
coordinates1[1,3] <- -83.8
coordinates1[1,2] <- 28.4

admixture_offset <- admixture1[admixture1$Site == "FL-7" ,]

admixture_offset <- standardise_data(admixture_offset, type = "admixture") |> transform_admix_data(data = _)
coordinates1 <- standardise_data(coordinates1, type = "coordinates")
admix_coords <- merge_coords_data(coordinates1, admixture_offset)

map2 + 
  geom_segment(data = coordinates2, 
               aes(x = Lon, y = Lat, xend = Lon1, yend = Lat1), 
               color = "black", size = 0.3, alpha = 1) +
  add_pie_charts(admix_coords,
                 admix_columns = 4:7,
                 lat_column = "lat",
                 lon_column = "lon",
                 pie_colours =  brewer.pal(4, "Set1"),
                 border = 0.5,
                 opacity = 0.9,
                 pie_size = 0.5
  )+
  theme(
    legend.title = element_blank(),
  )

# plot by snp
plot_snp <- function(locus) {
  snpfreqs <- filter(pop_freqs, loci == locus)
  snpfreqs <- rownames_to_column(snpfreqs, "site") 
  snpfreqs <- snpfreqs %>% arrange(site) %>%
    mutate(lat = arrange(coordinates, Site)$Lat) %>%
    mutate(lon = arrange(coordinates, Site)$Lon)
  par(fg = "black")
  maps::map("state", col = "grey85", fill = TRUE, border = FALSE, xlim=c(-100,-70), ylim=c(25,42))
  map.axes()
  for (i in 1:nrow(snpfreqs)){
    if (snpfreqs[i,2] + snpfreqs[i,3] != 0) {
      add.pie(z = c(snpfreqs[i,2], snpfreqs[i,3]), 
            x = snpfreqs$lon[i], 
            y = snpfreqs$lat[i], 
            radius = 0.5, col = c("blue","orange"), labels = "") 
    }
  }
}
snpfreqs <- filter(pop_freqs, loci == "Chr02_39132033")

cand_by_loading <- cand %>%
  arrange(loading)
plot_snp("Chr02_39132033")

snpfreqs <- pop_freqs %>%
  filter(loci == "Chr18_13842750")
snpfreqs <- rownames_to_column(snpfreqs, "site") 
snpfreqs <- snpfreqs %>% arrange(site) %>%
  mutate(lat = arrange(coordinates, Site)$Lat) %>%
  mutate(lon = arrange(coordinates, Site)$Lon)
par(fg = "black")
maps::map("state", col = "grey85", fill = TRUE, border = FALSE, xlim=c(-100,-70), ylim=c(25,42))
map.axes()
for (i in 1:nrow(snpfreqs)){
  if (snpfreqs[i,2] + snpfreqs[i,3] != 0) {
    add.pie(z = c(snpfreqs[i,2], snpfreqs[i,3]), 
            x = snpfreqs$lon[i], 
            y = snpfreqs$lat[i], 
            radius = 0.5, col = c("blue","orange"), labels = "") 
  }
}

for (candidate in range) {
  plot_snp(cand_by_loading$snp[candidate])
}

# Run mapmixture (chad version, by individual)
map2 <- mapmixture(
  admixture_df = admixture_by_ind,
  coords_df = coordinates,
  cluster_cols = c("#ff7f00","#1f78b4","#4c417a","#06592A","#ffac59","#59abe2","#8276b6","#0fe16a"),
  cluster_names = c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster1 weak","Cluster2 weak","Cluster3 weak","Cluster4 weak"),
  crs = 4326,
  boundary = c(xmin=-100, xmax=-72, ymin=26, ymax=38),
  #boundary = c(xmin=-90, xmax=-75, ymin=25, ymax=35),
  pie_size = 0.5,
  pie_border = 0.1,
  pie_border_col = "white",
  pie_opacity = 0.9,
  land_colour = "#d9d9d9",
  sea_colour = "#deebf7",
  basemap = rnaturalearthdata::states50,
  expand = TRUE,
  arrow = TRUE,
  arrow_size = 1.5,
  arrow_position = "bl",
  scalebar = TRUE,
  scalebar_size = 1.5,
  scalebar_position = "tl",
  plot_title = "Admixture Map",
  plot_title_size = 12,
  axis_title_size = 10,
  axis_text_size = 8
)
map2

# mapmixture + structure-style plot
map5 <- mapmixture(
  admixture_df = admixture1,
  coords_df = coordinates,
  cluster_cols = c("#ff7f00","#1f78b4","#998ec3","#06592A"),
  cluster_names = c("Group A","Group B","Group C","Group D"),
  crs = 4326,
  boundary = c(xmin=-100, xmax=-72, ymin=26, ymax=38),
  #boundary = c(xmin=-90, xmax=-75, ymin=25, ymax=35),
  basemap = rnaturalearthdata::states50,
  pie_size = 0.5,
  pie_border = 0.1,
  pie_border_col = "white",
  pie_opacity = 0.7,
  land_colour = "#d9d9d9",
  sea_colour = "#deebf7",
  expand = TRUE,
  arrow = TRUE,
  arrow_size = 1.5,
  arrow_position = "bl",
  scalebar = TRUE,
  scalebar_size = 1.5,
  scalebar_position = "tl",
  plot_title = "Admixture Map",
  plot_title_size = 12,
  axis_title_size = 10,
  axis_text_size = 8
) +
  # Adjust theme options
  theme(
    legend.position = "top",
    plot.margin = margin(l = 10, r = 10),
  ) +
  # Adjust the size of the legend keys
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1)))

# Traditional structure barplot
structure_barplot <- structure_plot(
  admixture_df = admixture1,
  type = "structure",
  cluster_cols = c("#ff7f00","#1f78b4","#998ec3","#06592A"),
  site_dividers = TRUE,
  divider_width = 0.4,
  site_order = coordinates$Site,
  labels = "site",
  flip_axis = FALSE,
  site_ticks_size = -0.05,
  site_labels_y = -0.35,
  site_labels_size = 2.2
)+
  # Adjust theme options
  theme(
    axis.title.y = element_text(size = 8, hjust = 1),
    axis.text.y = element_text(size = 5),
  )

# Arrange plots
grid.arrange(map5, structure_barplot, nrow = 2, heights = c(4,1))


#########################
### RDA ###
#########################

# import genetic data
setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1-5")

gen <- as(gen, "matrix")
sum(is.na(gen))
dim(gen)

gen[gen == 0] <- NA

# add pop cluster data
pops <- read.table("K4_str_clusters_Ivom.txt", header = TRUE)

clusters <- pops[,5:8]

pops[4] <- colnames(clusters)[apply(clusters,1,which.max)]
colnames(pops)[4] <- "clust"

gen.imp <- as.data.frame(gen)

sum(is.na(gen))

sum(is.na(gen.imp))
# remove markers with >20% missing data
na_pop <- apply(gen.imp[,-1], 2, function(x) sum(is.na(x)))
gen.imp <- gen.imp[,(which(na_pop<73)+1)]
print(paste0(length(gen.imp)-1," variants left after filtering"))

# impute missing genotypes based on median of subpopulation cluster
gen.imp <- add_column(gen.imp, pops[, 4], .before = 1)
colnames(gen.imp)[1] <- "pop"
gen_imp_by_group <- gen.imp %>% 
  dplyr::group_by(pop) %>%
  dplyr::mutate(across(dplyr::where(is.numeric), ~ coalesce(., median(., na.rm=TRUE)))) %>%
  ungroup()

gen_imp_by_group$pop <- NULL

# import climate data
bioclim <- worldclim_global(path="/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/worldclim", var="bio", res=2.5, version="2.1")
salinity <- rast("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/worldclim/HWSD2_RASTER/HWSD2.bil")
CMD <- rast("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/worldclim/cmd_1961-1990.tif")
sal <- terra::plot(salinity)
plot(bioclim[[14]])

# import coastline data
coast <- ne_coastline(scale=50, returnclass = "sf")

# can subset to only variables of interest
#bioclim <- bioclim[[c(3,19)]]

# read in sample data and clean rows with missing data
data <- read.csv("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1-5/Ivom1-5_wild_meta.csv", header=TRUE)


# subset out coordinates
coords <- data[, c(2,3)]   # coordinates

# convert to S4 object SpatVector with inherited coordinate system
# for bioclim and salinity
points <- vect(coords,
               geom=c("lon","lat"),
               crs = "EPSG:4326")

# get straight line distance to nearest coastline
dist2coast <- ggOceanMaps::dist2land(coords, shapefile=coast)

# extract climate data for points
values <- terra::extract(bioclim,points)
values2 <- terra::extract(salinity,points)
values3 <- terra::extract(CMD, points)
# add other data back in
clim <- cbind.data.frame(data,values,values2,values3,dist2coast[3])

env <- clim[, -c(10,30,32)]

env$id <- as.character(env$id)
identical(rownames(gen), env[,1]) 
# find filtered individuals and remove
# which(env$id %in% c("FL-MZ-4","FL-MZ-5","AR-CT-21","BG-2"))
# env <- env %>% slice(-c(7,18,409,411))
# identical(rownames(gen), env[,1]) 


# transform coast data
hist(env$ldist)
env <- env %>%
  mutate(logCoast = log2(ldist))
hist(env$logCoast)

# add first three PCs from PCA on genetic data using intergenic regions
env <- env %>%
  mutate(labelcheck = pcs$.rownames,
         PC1 = pcs$.fittedPC1,
         PC2 = pcs$.fittedPC2,
         PC3 = pcs$.fittedPC3,
         )
identical(env$id, env$labelcheck) 
env$labelcheck <- NULL

#####
## Various Model Building and Exploring Autocorrelation Between Predictors
#####

# variable selection using forward model building procedure
RDA0 <- rda(formula = gen_imp_by_group ~ 1, data = env, scale = T) 
words <- colnames(env)
words1 <- strsplit(words, "\\s+")
variables <- paste0(words1, collapse="+")

# Stepwise model building with ordiR2step
var_selection <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)

#examine correlation between predictors
pairs.panels(env[,9:29], scale=TRUE)

######
## select environmental variables for RDA and check colinearity
######

#pred <- subset(env, select=-c(sex,cmt,wc2.1_2.5m_bio_1,wc2.1_2.5m_bio_4,wc2.1_2.5m_bio_7,wc2.1_2.5m_bio_10,wc2.1_2.5m_bio_11,wc2.1_2.5m_bio_15,wc2.1_2.5m_bio_16,wc2.1_2.5m_bio_17,wc2.1_2.5m_bio_18,wc2.1_2.5m_bio_19))
pred <- subset(env, select=-c(wc2.1_2.5m_bio_1,wc2.1_2.5m_bio_2,wc2.1_2.5m_bio_3,wc2.1_2.5m_bio_4,wc2.1_2.5m_bio_7,wc2.1_2.5m_bio_8,wc2.1_2.5m_bio_9,wc2.1_2.5m_bio_10,wc2.1_2.5m_bio_11,wc2.1_2.5m_bio_12,wc2.1_2.5m_bio_13,wc2.1_2.5m_bio_15,wc2.1_2.5m_bio_16,wc2.1_2.5m_bio_17,wc2.1_2.5m_bio_18,wc2.1_2.5m_bio_19))

pred <- pred %>%
  relocate(ele, .after = superpop)

pairs.panels(pred[,9:15], scale=TRUE)

## set column names
#colnames(pred) <- c("ID","lat","lon","ele","MDR","IT","MAXWM","MINCM","MTWQ","MTDQ","AP","PWM","PDM")
colnames(pred) <- c("id","lat","lon","sex","site","state","ecotype","population","ele","MAXWM","MINCM","PDM","SAL","CMD","COAST","PC1","PC2","PC3")

#############
## the RDA
#############
# Exploring variance partitioning - full model
RDAfull <- rda(formula = gen ~ PC1 + PC2 + PC3 + lat + lon + MAXWM + MINCM + PDM + SAL + COAST + CMD, data = pred, scale = T)
# only genetic structure
pRDAstruct <- rda(formula = gen ~ PC1 + PC2 + PC3 + Condition(lat + lon + MAXWM + MINCM + PDM + SAL + COAST + CMD), data = pred, scale = T)
# only geography
pRDAgeo <- rda(formula = gen ~ lat + lon + Condition(PC1 + PC2 + PC3 + MAXWM + MINCM + PDM + SAL + COAST + CMD), data = pred, scale = T)
# climate-only model
Ivom.rda <- rda(formula = gen ~ ele + MAXWM + MINCM + PDM + SAL + COAST + CMD + Condition(PC1 + PC2 + PC3), data = pred, scale = T)
RsquareAdj(RDAfull)
RsquareAdj(pRDAstruct)
RsquareAdj(pRDAgeo)
RsquareAdj(Ivom.rda)
RsquaredSummary <- list(RsquareAdj(RDAfull),RsquareAdj(pRDAstruct),RsquareAdj(pRDAgeo),RsquareAdj(Ivom.rda))
v1 <- c("Full Model","Genetic Structure Only","Georgraphy Only","Climate Only")
v2 <- as_vector(c(RsquaredSummary[[1]][[1]], RsquaredSummary[[2]][[1]],
                RsquaredSummary[[3]][[1]], RsquaredSummary[[4]][[1]]))
v3 <- as_vector(c(RsquaredSummary[[1]][[2]], RsquaredSummary[[2]][[2]],
                  RsquaredSummary[[3]][[2]], RsquaredSummary[[4]][[2]]))
RsquaredSummary <- cbind(v1,v2,v3)
colnames(RsquaredSummary) <- c("Model", "R_squared", "R_squared_adj")
as_tibble(RsquaredSummary)
write.csv(RsquaredSummary, "RDA_Rsq_values.csv", row.names = FALSE)

summary(eigenvals(Ivom.rda, model = "constrained"))
signif.full <- anova.cca(Ivom.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full
signif.axis <- anova.cca(Ivom.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
vif.cca(Ivom.rda)

screeplot(Ivom.rda) # determine how many axes to use (K)
Ivom.rdadapt <- rdadapt(RDA = Ivom.rda, K = 3)
# P-values threshold after Bonferroni correction
thres <- 0.01/length(Ivom.rdadapt$p.values)
# Identifying the loci that are below the p-value threshold
outliers <- data.frame(Loci = colnames(gen)[which(Ivom.rdadapt$p.values<thres)], p.value = Ivom.rdadapt$p.values[which(Ivom.rdadapt$p.values<thres)], contig = unlist(lapply(strsplit(colnames(gen)[which(Ivom.rdadapt$p.values<thres)], split = "_"), function(x) x[1])))
length(outliers)
plot(Ivom.rda, scaling=3)
plot(Ivom.rda, choices = c(1, 3), scaling=3)

levels(pred$population) <- c("atlantic","gulf","florida")
#levels(env$pop) <- c("AR","NC","VA","FL","GA","SC","AL","MS","TX","LA","BW","MC")
eco <- env$population
#eco <- env$pop
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#9CCB86","#F2ACCA","#A6761D","#AF58BA","#009392","#7C1D6F","#B10026","#06592A") # 6 nice colors for our ecotypes
#levels(env$gen_clust) <- c("c1","c2","c3","c4")
#gen_clusts <- env$gen_clust
bg <- c("#1f78b4","#4c417a","#06592A")

# axes 1 & 2
plot(Ivom.rda, type="n", scaling=3)
points(Ivom.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(Ivom.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg) # the plants
text(Ivom.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(Ivom.rda, type="n", scaling=3, choices=c(1,3))
points(Ivom.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(Ivom.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg, choices=c(1,3))
text(Ivom.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

load.rda <- scores(Ivom.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3.5) # 
cand2 <- outliers(load.rda[,2],3.5) # 
cand3 <- outliers(load.rda[,3],3.5) # 

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=7)  # 7 columns for 7 predictors
colnames(foo) <- c("ele","MAXWM","MINCM","PDM","SAL","CMD","COAST")

## calculate correlation of each candidate SNP to env variables
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen[,nam]
  foo[i,] <- apply(pred[, c(9,10,11,12,13,14,15)],2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

## check for duplicates
length(cand$snp[duplicated(cand$snp)])
cand <- cand[!duplicated(cand$snp),]

#find which predictors each candidate SNP is most associated with
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,11] <- names(which.max(abs(bar[4:10]))) # gives the variable
  cand[i,12] <- max(abs(bar[4:10]))              # gives the correlation
}

colnames(cand)[11] <- "predictor"
colnames(cand)[12] <- "correlation"

table(cand$predictor) 

#plot SNPs 
sel <- cand$snp
env <- cand$predictor
#env[env=="ele"] <- '#1f78b4'
env[env=="ele"] <- '#DC69B9'
env[env=="MAXWM"] <- '#a6cee3'
env[env=="MINCM"] <- '#6a3d9a'
env[env=="PDM"] <- '#e31a1c'
env[env=="SAL"] <- '#33a02c'
env[env=="CMD"] <- 'darkorange'
env[env=="COAST"] <- '#ffff33'

# color by predictor:
col.pred <- rownames(Ivom.rda$CCA$v) # pull all the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("Chr",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c',"darkorange",'#DC69B9')

# axes 1 & 2
plot(Ivom.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
# plot all snps with nonsig greyed out
#points(Ivom.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
# plot only sig snps
points(Ivom.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(Ivom.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("topright", legend=c("ele","MAXWM","MINCM","PDM","SAL","CMD","COAST"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(Ivom.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
#points(Ivom.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(Ivom.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(Ivom.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topright", legend=c("MAXWM","MINCM","PDM","SAL","CMD","COAST"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 2 & 3
plot(Ivom.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3))
#points(Ivom.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(Ivom.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(Ivom.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topright", legend=c("MAXWM","MINCM","PDM","SAL","CMD","COAST"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

write.csv(cand, "./results_expanded_data/rda_cand_snps.csv")
###########################
## LD decay plotting script 
###########################

rm(list = ls())
install.packages("tidyverse")
library(tidyverse)

# set path
my_bins <- "./Ivom_chr1.ld_decay_bins"

# read in data
ld_bins <- read_table(my_bins)

# plot LD decay
ggplot(ld_bins, aes(distance, avg_R2)) + geom_line() +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2)) +
  coord_cartesian(xlim = c(0,250000))

##########

popquery <- read_delim("popquery", delim=">", col_names=F)
popcommand <- popquery %>% mutate(paste0("--weir-fst-pop", popquery$X2, " \\"))
colnames(popcommand) <- c("X1","X2","X3")
command <- popcommand %>% dplyr::select(X3)
write_tsv(command, "commands") 

#import Fst data
fst <- read_tsv("Ivom384allsites.weir.fst")

# capdown the headers
names(fst) <- tolower(names(fst))
colnames(fst)[3] <- "fst"

ggplot(fst, aes(pos, fst)) + geom_point()

# identify the 95% and 99% percentile
quantile(fst$fst, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(fst$fst, 0.975, na.rm = T)
# make an outlier column in the data.frame
fst <- fst %>% mutate(outlier = ifelse(fst > my_threshold, "outlier", "background"))
fst %>% group_by(outlier) %>% tally()
ggplot(fst, aes(pos, fst, colour = outlier)) + geom_point()

ggplot(fst,aes(fst)) + geom_histogram()

## utility script for sorting sample names ##

id_sex <- read.table("id-sex.txt", sep="\t")
id_sex <- id_sex[-1,]
sample_order <- read.table("sample_order.txt")
id_sex$V1 <- factor(id_sex$V1, levels = sample_order$V1)
id_sex <- id_sex[order(id_sex$V1),]
write.csv(id_sex, "id_sex_sorted.csv")

