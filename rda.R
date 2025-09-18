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
library(mapmixture)
library(geodata)
library(raster)
library(sp)
library(terra)
library(psych)
library(vegan)
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
## plot structure results onto map

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
# Run mapmixture
map1 <- mapmixture(admixture1, coordinates, crs = 3035)
map1

# Run mapmixture (chad version, by mean membership per pop)
map2 <- mapmixture(
  admixture_df = admixture1,
  coords_df = coordinates,
  cluster_cols = c("#ff7f00","#1f78b4","#4c417a","#06592A"),
  cluster_names = c("Cluster1","Cluster2","Cluster3","Cluster4"),
  crs = 4326,
  boundary = c(xmin=-100, xmax=-72, ymin=26, ymax=38),
  #boundary = c(xmin=-90, xmax=-75, ymin=25, ymax=35),
  pie_size = 0.5,
  pie_border = 0.1,
  pie_border_col = "black",
  pie_opacity = 0.75,
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
                 pie_colours = c("#ff7f00","#1f78b4","#4c417a","#06592A"),
                 border = 0.1,
                 opacity = 1,
                 pie_size = 0.5
  )+
  theme(
    legend.title = element_blank(),
  )



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

### RDA ###

#import genetic data

gen <- as(Ivom384, "matrix")
sum(is.na(gen))
dim(gen)

gen[gen == 0] <- NA

sum(is.na(gen.imp)) # No NAs

# add pop cluster data
pops <- read.table("K4_str_clusters_Ivom", header = TRUE)

clusters <- pops[,5:8]

pops[4] <- colnames(clusters)[apply(clusters,1,which.max)]
colnames(pops)[4] <- "clust"

gen.imp <- as.data.frame(gen)

gen.imp <- add_column(gen.imp, pops[, 4], .before = 1)
colnames(gen.imp)[1] <- "pop"

# remove markers with >20% missing data
na_pop <- apply(gen.imp[,-1], 2, function(x) sum(is.na(x)))
gen.imp <- gen.imp[,(which(na_pop<73)+1)]
print(paste0(length(gen.imp)-1," variants left after filtering"))

# impute missing genotypes based on median of subpopulation cluster
gen_imp_by_group <- gen.imp %>% 
  group_by(pop) %>%
  mutate(across(where(is.numeric), ~ coalesce(., median(., na.rm=TRUE)))) %>%
  ungroup()

gen_imp_by_group$pop <- NULL

# import climate data
bioclim <- worldclim_global(path="/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/worldclim", var="bio", res=2.5, version="2.1")
salinity <- rast("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/worldclim/HWSD2_RASTER/HWSD2.bil")
sal <- terra::plot(salinity)

# import coastline data
coast <- ne_coastline(scale=50, returnclass = "sf")

# can subset to only variables of interest
#bioclim <- bioclim[[c(3,19)]]

# read in sample data and clean rows with missing data
data <- read.csv("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/worldclim/yaupon-diversity-panel-ordered.csv", header=TRUE)


# subset out coordinates
coords <- data[, c(3,4)]   # coordinates

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

# add other data back in
clim <- cbind.data.frame(data,values,values2,dist2coast[3])

env <- clim[, -c(9,29)]

env$ID <- as.character(env$ID)
identical(rownames(gen), env[,1]) 

# variable selection using forward model building procedure
RDA0 <- rda(formula = gen_imp_by_group ~ 1, data = env, scale = T) 
words <- colnames(env)
words1 <- strsplit(words, "\\s+")
variables <- paste0(words1, collapse="+")

RDAfull <- rda(gen_imp_by_group ~ lat+lon+ele+wc2.1_2.5m_bio_1+wc2.1_2.5m_bio_2+wc2.1_2.5m_bio_3+wc2.1_2.5m_bio_4+wc2.1_2.5m_bio_5+wc2.1_2.5m_bio_6+wc2.1_2.5m_bio_7+wc2.1_2.5m_bio_8+wc2.1_2.5m_bio_9+wc2.1_2.5m_bio_10+wc2.1_2.5m_bio_11+wc2.1_2.5m_bio_12+wc2.1_2.5m_bio_13+wc2.1_2.5m_bio_14+wc2.1_2.5m_bio_15+wc2.1_2.5m_bio_16+wc2.1_2.5m_bio_17+wc2.1_2.5m_bio_18+wc2.1_2.5m_bio_19+HWSD2+ldist, data=env, scale=T)

# Stepwise model building with ordiR2step
var_selection <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)

#examine correlation between predictors
pairs.panels(env[,9:29], scale=TRUE)

#pred <- subset(env, select=-c(sex,cmt,wc2.1_2.5m_bio_1,wc2.1_2.5m_bio_4,wc2.1_2.5m_bio_7,wc2.1_2.5m_bio_10,wc2.1_2.5m_bio_11,wc2.1_2.5m_bio_15,wc2.1_2.5m_bio_16,wc2.1_2.5m_bio_17,wc2.1_2.5m_bio_18,wc2.1_2.5m_bio_19))
pred <- subset(env, select=-c(cmt,wc2.1_2.5m_bio_1,wc2.1_2.5m_bio_2,wc2.1_2.5m_bio_3,wc2.1_2.5m_bio_4,wc2.1_2.5m_bio_7,wc2.1_2.5m_bio_8,wc2.1_2.5m_bio_9,wc2.1_2.5m_bio_10,wc2.1_2.5m_bio_11,wc2.1_2.5m_bio_12,wc2.1_2.5m_bio_13,wc2.1_2.5m_bio_15,wc2.1_2.5m_bio_16,wc2.1_2.5m_bio_17,wc2.1_2.5m_bio_18,wc2.1_2.5m_bio_19))


pairs.panels(pred[,4:7], scale=TRUE)

#colnames(pred) <- c("ID","lat","lon","ele","MDR","IT","MAXWM","MINCM","MTWQ","MTDQ","AP","PWM","PDM")
colnames(pred) <- c("ID","pop","lat","lon","ele","sex","gen_clust","MAXWM","MINCM","PDM","SAL","COAST")

#Ivom384.rda <- rda(formula = gen.imp ~ ele + MDR + IT + MAXWM + MINCM + MTWQ + MTDQ + AP + PWM + PDM, data = pred, scale = T)
pRDA_full <- rda(formula = gen_imp_by_group ~ ele + lat + lon + MAXWM + MINCM + PDM + SAL + COAST, data = pred, scale = T)
RsquareAdj(pRDA_full)

Ivom384.rda <- rda(formula = gen_imp_by_group ~ ele + MAXWM + MINCM + PDM + SAL + COAST, data = pred, scale = T)

RsquareAdj(Ivom384.rda)
summary(eigenvals(Ivom384.rda, model = "constrained"))
signif.full <- anova.cca(Ivom384.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full
signif.axis <- anova.cca(wolf.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
vif.cca(Ivom384.rda)

plot(Ivom384.rda, scaling=3)
plot(Ivom384.rda, choices = c(1, 3), scaling=3)

levels(env$pop) <- c("AR","NC","VA","FL","GA","SC","AL","MS","TX","LA","BW","MC")
eco <- env$pop
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#9CCB86","#F2ACCA","#A6761D","#AF58BA","#009392","#7C1D6F","#B10026","#06592A") # 6 nice colors for our ecotypes

levels(env$gen_clust) <- c("c1","c2","c3","c4")
gen_clusts <- env$gen_clust
bg <- c("#ff7f00","#1f78b4","#4c417a","#06592A")

# axes 1 & 2
plot(Ivom384.rda, type="n", scaling=3)
points(Ivom384.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(Ivom384.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg) # the plants
text(Ivom384.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(Ivom384.rda, type="n", scaling=3, choices=c(1,3))
points(Ivom384.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(Ivom384.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg, choices=c(1,3))
text(Ivom384.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

load.rda <- scores(Ivom384.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3.5) # 38
cand2 <- outliers(load.rda[,2],3.5) # 69
cand3 <- outliers(load.rda[,3],3.5) # 34

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=6)  # 6 columns for 6 predictors
colnames(foo) <- c("ele","MAXWM","MINCM","PDM","SAL","COAST")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen_imp_by_group[,nam]
  foo[i,] <- apply(pred[, c(5,8,9,10,11,12)],2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$snp[duplicated(cand$snp)])
cand <- cand[!duplicated(cand$snp),]

#find which predictors each candidate SNP is most associated with
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,10] <- names(which.max(abs(bar[4:9]))) # gives the variable
  cand[i,11] <- max(abs(bar[4:9]))              # gives the correlation
}

colnames(cand)[10] <- "predictor"
colnames(cand)[11] <- "correlation"

table(cand$predictor) 

#plot SNPs 
sel <- cand$snp
env <- cand$predictor
env[env=="ele"] <- '#1f78b4'
env[env=="MAXWM"] <- '#a6cee3'
env[env=="MINCM"] <- '#6a3d9a'
env[env=="PDM"] <- '#e31a1c'
env[env=="SAL"] <- '#33a02c'
env[env=="COAST"] <- '#ffff33'


# color by predictor:
col.pred <- rownames(Ivom384.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("h1tg",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33')

# axes 1 & 2
plot(Ivom384.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(Ivom384.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(Ivom384.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(Ivom384.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("ele","MAXWM","MINCM","PDM","SAL","COAST"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(Ivom384.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(Ivom384.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(Ivom384.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(Ivom384.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c("ele","MAXWM","MINCM","PDM","SAL","COAST"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
