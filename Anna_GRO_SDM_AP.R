setwd("~/Desktop/Granite Outcrops/GRO_SDM")
library(raster)
library(sp)
library(maptools)
library(rgdal)
library(dismo)
library(sf)
library(rJava)
library(envirem)

# SDM for Gratiola amphiantha
# known presences for Gratiola amphiantha
ap.occ <- read.csv("GRO_AllKnownSitesEver_AP.csv")
ap.occ <- ap.occ[-c(1,2,5:12)]
plot(ap.occ$Long,ap.occ$Lat,xlab='Longitude',ylab='Latitude')

# establish spatial extent of the study
max.lat <- ceiling(35.10030)
min.lat <- floor(32.68413)
max.lon <- ceiling(-80.52640)
min.lon <- floor(-85.69899)

# create a map of the occurrences
data(wrld_simpl)

plot(wrld_simpl,xlim=c(min.lon, max.lon),
     ylim=c(min.lat, max.lat),axes=T, col='grey95')
points(x=ap.occ$Long, y=ap.occ$Lat, col='olivedrab',
       pch=20, cex=0.75)

### load environmental data from worldclim

# first create a folder in your working directory called "data"
# you'll download the environmental data into this folder

# note: you need to set download=TRUE the first time you download these files
# but you can change to download=F for any subsequent time you re-fun the script
bio_curr <- getData("worldclim", var="bio", res=2.5, download=TRUE, path="data/")
bio_fut <- getData("CMIP5", var='bio', res=2.5, 
                   rcp=45, year=50, model='NO', download=TRUE, path="data/")
bio_past <- getData("envirem", var="bio", res-2.5, download = TRUE, path = "data/")

# crop climate data to the spatial extent of the study
study <- extent(min.lon, max.lon, min.lat, max.lat)
bio_curr2 <- crop(bio_curr, study)
bio_fut2 <- crop(bio_fut, study)
# give the climate variables the same names
names(bio_fut2)<-names(bio_curr2)
# visualize current and future environmental data
plot(bio_curr2)
plot(bio_fut2)

### generate background points (pseudo-absence data)
# use the bioclim files to get spatial extent and sampling resolution
bio.files <- list.files(path="data/wc2-5", pattern="*.bil$",full.names=T)
bio.files
mask <- raster(bio.files[1]) # use the first file to define spatial extent and scale

bg.points <- randomPoints(mask=mask, n=nrow(ap.occ), excludep=TRUE, ext=study, extf = 0.99)
head(bg.points)

# visualize background points
plot(wrld_simpl,xlim=c(min.lon, max.lon),
     ylim=c(min.lat, max.lat),axes=T, col='grey95')
points(x=ap.occ$Long, y=ap.occ$Lat, col='olivedrab',
       pch=20, cex=0.75)
points(bg.points, col='grey30',
       pch=1, cex=0.75)

# now combine into one dataframe
bg.points <- as.data.frame(bg.points)
colnames(bg.points) <- c('Long', 'Lat')
ap.occ$presence <- 1
bg.points$presence <- 0
all_points <- rbind(ap.occ, bg.points)
head(all_points)

### split into testing and training datasets for cross-validation

# we'll arbitrarily set the "testing group" to 1
testing.group <- 1

# create random groups (traditional)
random.group <- kfold(x=all_points, k=10)
table(random.group) # number of observations per group
random.train <- all_points[random.group != testing.group,]
random.test <- all_points[random.group == testing.group,]

############# Fit some distribution models!
# Here we're using Maxent
# fit the model
mod.me <- maxent(x=bio_curr2, 
                 p=random.train[random.train$presence==1,c('Long','Lat')],
                 a=random.train[random.train$presence==0,c('Long','Lat')])
# look at the predictions based on current climate
mod.me.predict <- predict(mod.me, bio_curr2)
plot(mod.me.predict, main='Current climate')

# look at the importance of each climate variable
plot(mod.me) # it's basically all bio6!, some bio8,4,14,10,18

# look at the partial response curves to each climate variable
# (holding others constant at their mean)
response(mod.me)

# evaluate how well the model does with the testing dataset
mod.me.eval <- evaluate(mod.me, x=bio_curr2,
                        p=random.test[random.test$presence==1,c('Long','Lat')],
                        a=random.test[random.test$presence==0,c('Long','Lat')])
mod.me.eval@auc # auc > 0.7 considered ok, this is 0.6
plot(mod.me.eval,'ROC')

# look at the predictions based on future climate
mod.me.forecast <- predict(mod.me, x=bio_fut2, ext=study)
plot(mod.me.forecast, main='RCP 4.5')

# convert continuous suitability scores into binary presence/absence maps
thresh <- threshold(mod.me.eval, stat='spec_sens')
plot(mod.me.predict>thresh, main='Current climate')
plot(mod.me.forecast>thresh, main='RCP 4.5')

# need to extract prediction values from observed EOs
# import individual pool lat longs, export SDM values from these coordinates

############ SDM from Fletcher & Fortin ############
library(raster)
GA <- raster()

