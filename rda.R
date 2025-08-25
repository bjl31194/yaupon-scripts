install.packages("raster")
install.packages("sp")
install.packages("geodata")
install.packages("psych")
library(geodata)
library(raster)
library(sp)
library(terra)
library(psych)
library(vegan)

#import genetic data

gen <- as(Ivom384, "matrix")
sum(is.na(gen))
dim(gen)

gen[gen == 0] <- NA
# ONLY FOR TESTING, impute missing genotypes (replace NAs with most common genotype at each locus)
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs

# import climate data
bioclim <- worldclim_global(path="/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/worldclim", var="bio", res=2.5, version="2.1")

# can subset to only variables of interest
#bioclim <- bioclim[[c(3,19)]]

# read in sample data and clean rows with missing data
data <- read.csv("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/worldclim/yaupon-diversity-panel-ordered.csv", header=TRUE)


# subset out coordinates
coords <- data[, c(3,4)]   # coordinates

# convert to S4 object SpatVector with inherited coordinate system
points <- vect(coords,
               geom=c("lon","lat"),
               crs = bioclim@pntr)

# extract climate data for points
values <- terra::extract(bioclim,points)

# add other data back in
clim <- cbind.data.frame(data,values)

env <- clim[, -8]

env$ID <- as.character(env$ID)
identical(rownames(gen), env[,1]) 

pairs.panels(env[,7:25], scale=TRUE)

#pred <- subset(env, select=-c(sex,cmt,wc2.1_2.5m_bio_1,wc2.1_2.5m_bio_4,wc2.1_2.5m_bio_7,wc2.1_2.5m_bio_10,wc2.1_2.5m_bio_11,wc2.1_2.5m_bio_15,wc2.1_2.5m_bio_16,wc2.1_2.5m_bio_17,wc2.1_2.5m_bio_18,wc2.1_2.5m_bio_19))
pred <- subset(env, select=-c(cmt,wc2.1_2.5m_bio_1,wc2.1_2.5m_bio_2,wc2.1_2.5m_bio_3,wc2.1_2.5m_bio_4,wc2.1_2.5m_bio_7,wc2.1_2.5m_bio_8,wc2.1_2.5m_bio_9,wc2.1_2.5m_bio_10,wc2.1_2.5m_bio_11,wc2.1_2.5m_bio_12,wc2.1_2.5m_bio_13,wc2.1_2.5m_bio_15,wc2.1_2.5m_bio_16,wc2.1_2.5m_bio_17,wc2.1_2.5m_bio_18,wc2.1_2.5m_bio_19))


pairs.panels(pred[,4:7], scale=TRUE)

#colnames(pred) <- c("ID","lat","lon","ele","MDR","IT","MAXWM","MINCM","MTWQ","MTDQ","AP","PWM","PDM")
colnames(pred) <- c("ID","pop","lat","lon","ele","sex","MAXWM","MINCM","PDM")


#Ivom384.rda <- rda(formula = gen.imp ~ ele + MDR + IT + MAXWM + MINCM + MTWQ + MTDQ + AP + PWM + PDM, data = pred, scale = T)
Ivom384.rda <- rda(formula = gen.imp ~ ele + MAXWM + MINCM + PDM, data = pred, scale = T)

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

