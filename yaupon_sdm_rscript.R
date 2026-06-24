setwd()
library(raster)
library(terra)
library(sp)
install.packages("maptools")
library(maptools)
library(rgdal)
library(dismo)
library(sf)
install.packages("rJava")
library(rJava)
library(envirem)
install.packages("spocc")
library(spocc)
install.packages("PresenceAbsence")
library(PresenceAbsence)
install.packages("DAAG")
library(DAAG)
library(spData)

# read in occurrence data
iv.occ <- read_csv("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/sdm/occurrences/gbif_occ_RAW.csv")
plot(iv.occ$x,iv.occ$y,xlab='decimalLongitude',ylab='decimalLatitude')

# establish spatial extent of the study

max.lat <- ceiling(42)
min.lat <- floor(24)
max.lon <- ceiling(-75)
min.lon <- floor(-105)
study <- extent(min.lon, max.lon, min.lat, max.lat)

# create a continental map for cropping
usa_states <- st_as_sf(spData::us_states)
continental_us <- usa_states[usa_states$NAME != "Alaska" & usa_states$NAME != "Hawaii", ]

map <- ne_countries(
  continent = "North America",
  type = "map_units", scale = "medium"
)


# Ensure coordinate reference systems match
continental_us <- st_transform(continental_us, crs(bio_past_ext))

### load environmental data from worldclim

# note: you need to set download=TRUE the first time you download these files
# but you can change to download=F for any subsequent time you re-fun the script
bio_pres <- worldclim_global(path="/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/worldclim/climate/wc2.1_2.5m", download = FALSE, var="bio", res=2.5, version="2.1")

terra::plot(bio_past_ext[[1]])
files <- list.files(path = "/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/worldclim/chelsa", pattern = 'tif',full.names = TRUE)
bio_past <- raster::stack(files)

# crop climate data to the spatial extent of the study

bio_pres_ext <- crop(bio_pres, study)
bio_past_ext <- crop(bio_past, study)
bio_past_ext <- crop(bio_past, bio_pres_ext)


terra::plot(bio_past_ext)
bio_past_ext <- na.omit(bio_past_ext)
bio_past_ext_summary <- summary(bio_past_ext)
bio_pres_ext_summary <- summary(bio_pres_ext)
View(bio_pres_ext_summary)
View(bio_past_ext_summary)
# give the climate variables the same names and ensure they have the same crs
names(bio_past_ext) <- names(bio_pres_ext)
crs(bio_past_ext) == crs(bio_pres_ext)
# convert units on LGM data to Worldclim (present) units



bio_past_ext[[1]] <- bio_past_ext[[1]] - 273.15
bio_past_ext[[5]] <- bio_past_ext[[5]] - 273.15
bio_past_ext[[6]] <- bio_past_ext[[6]] - 273.15
bio_past_ext[[7]] <- bio_past_ext[[7]] - 273.15
bio_past_ext[[8]] <- bio_past_ext[[8]] - 273.15
bio_past_ext[[9]] <- bio_past_ext[[9]] - 273.15
bio_past_ext[[10]] <- bio_past_ext[[10]] - 273.15
bio_past_ext[[11]] <- bio_past_ext[[11]] - 273.15
  
# visualize current and future environmental data

### generate background points (pseudo-absence data)
# use the bioclim files to get spatial extent and sampling resolution
bio.files <- list.files(path="/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/worldclim/climate/wc2.1_2.5m", pattern="*.tif",full.names=T)
bio.files
mask <- raster(bio.files[1]) # use the first file to define spatial extent and scale

bg.points <- randomPoints(mask=mask, n=nrow(iv.occ), excludep=TRUE, ext=study, extf = 0.99)
bg.points <- as.data.frame(bg.points)
head(bg.points)

# visualize background points
iv.occ <- st_as_sf(iv.occ, coords=c("decimalLongitude", "decimalLatitude"), crs = 4326)
bg.points <- st_as_sf(bg.points, coords=c("x", "y"), crs = 4326)
ggplot() +
  geom_sf(data = map) +
  geom_sf(data = iv.occ$geometry, color = "olivedrab") +
  geom_sf(data = bg.points$geometry, color = "black") +
  coord_sf(
    xlim = c(-105, -75),
    ylim = c(20, 42)
  ) +
  theme_minimal()

ivo_coords <- st_coordinates(iv.occ)
ivo <- cbind(ivo_coords, iv.occ)
ivo$PA <- 1
ivo <- ivo[c(1,2,51,52)]

abs_coords <- st_coordinates(bg.points)
ivo_abs <- cbind(abs_coords, bg.points)
ivo_abs$PA <- 0

ivo_PA <- rbind(ivo, ivo_abs) # combine 
# create SpatVector from PA data and extract env variables
spec.pres.sv <- as(ivo_PA, "SpatVector") # presence must be spatvector
spec.pres.extract <- terra::extract(bio_pres, spec.pres.sv) # perform extraction of preds to presence

spec.trained.data <- cbind(spec.pres.extract, ivo_PA) # combine true pres, fnetids, and extracted preds
head(spec.trained.data) # check
spec.trained.data.max <- na.exclude(spec.trained.data)

spec.trained.data <- na.omit(spec.trained.data)


# GLM Model Formula
mod.form <- function(dat, r.col, p.col){
  n.col <- ncol(dat) # Number of columns in the dataframe
  resp <- colnames(dat[r.col]) # assign response a column name
  pred <- colnames(dat[c(p.col:n.col)]) # assign preds column names
  mod.formula <- as.formula(paste(resp,
                 "~", paste(pred, collapse = "+"))) # formula
}

# Basic GLM with link = Binomial, dataframe = spec.trained.data
mod1.LR <- glm(as.factor(PA) ~ wc2.1_2.5m_bio_1 + wc2.1_2.5m_bio_12 + wc2.1_2.5m_bio_3 +
                 wc2.1_2.5m_bio_4 + wc2.1_2.5m_bio_5 + wc2.1_2.5m_bio_6 + wc2.1_2.5m_bio_7 + wc2.1_2.5m_bio_8 +
                 wc2.1_2.5m_bio_9 + wc2.1_2.5m_bio_10 + wc2.1_2.5m_bio_11 + wc2.1_2.5m_bio_12 + wc2.1_2.5m_bio_13 +
                 wc2.1_2.5m_bio_14 + wc2.1_2.5m_bio_15 + wc2.1_2.5m_bio_16 + wc2.1_2.5m_bio_17 + wc2.1_2.5m_bio_18 +
                 wc2.1_2.5m_bio_19, 
               family = binomial, data = spec.trained.data)

summary(mod1.LR)

# model 1 fit
mod1.fit <- 100 * (1 - mod1.LR$deviance/mod1.LR$null.deviance) # model fit
mod1.fit  # examine

mod1.pred <- predict(mod1.LR, type = "response") # model prediction
head(mod1.pred) # examine prediction

# LR model 2; backwards stepwise variable reduction
mod2.LR <- step(mod1.LR, trace = F)
# model 2 fit
100 * (1 - mod2.LR$deviance/mod2.LR$null.deviance)

# model 2 prediction
mod2.pred <- predict(mod2.LR, type = "response")
head(mod2.pred) # model 2
# model 2 summary
summary(mod2.LR)

# add var to keep track of model
modl <- "mod2.LR"
dat2 <- cbind(modl, spec.trained.data[23], mod2.pred) # build prediction df with mod2 preds

# Create confusion matrix
mod.cut <- optimal.thresholds(dat2, opt.methods = c("MaxKappa"))
mod2.cfmat <- table(dat2[[2]],
                    factor(as.numeric(dat2$mod2.pred >= mod.cut$mod2.pred)))
mod2.cfmat # examine confusion matrix

# Model accuracies 
mod2.acc <- presence.absence.accuracy(dat2, 
                                      threshold = mod.cut$mod2.pred,
                                      st.dev = F)
tss <- mod2.acc$sensitivity + mod2.acc$specificity - 1 
mod2.acc <- cbind(mod2.acc[1:7], tss) # bind metrics
mod2.acc[c(1, 4:5, 7:8)] # examine accuracies

auc.roc.plot(dat2, color = T) # visualize accuracy

# Cross validation accuracies
mod2.cv10 <- CVbinary(mod2.LR, nfolds = 5, print.details = F) # crossval
ls(mod2.cv10)

mod2.cv10.1 <- mod2.cv10$cvhat

dat2 <- cbind(dat2, mod2.cv10.1)
head(dat2)

# model CVbinary
mod2.cfmatCV <- table(dat2[[2]], 
                      factor(as.numeric(dat2$mod2.cv10.1 >= mod.cut$mod2.pred)))
# examine both mod confusion matrices
mod2.cfmatCV; mod2.cfmat

# calculate accuracies with std.dev = F
mod2.accB <- presence.absence.accuracy(dat2, 
                                       threshold = mod.cut$mod2.pred, 
                                       st.dev = F)
tss <- mod2.accB$sensitivity + mod2.accB$specificity - 1 # code TSS metric
mod2.accB <- cbind(mod2.accB[1:7], tss) # bind all metrics
mod2.accB[c(1, 4:5, 7:8)] # examine accuracies

auc.roc.plot(dat2, color = T)

# Cross validation accuracies
mod2.cv10 <- CVbinary(mod2.LR, nfolds = 5, print.details = F) # crossval
ls(mod2.cv10)
mod2.cv10.1 <- mod2.cv10$cvhat
dat2 <- cbind(dat2, mod2.cv10.1)
head(dat2)

# model CVbinary
mod2.cfmatCV <- table(dat2[[2]], 
                      factor(as.numeric(dat2$mod2.cv10.1 >= mod.cut$mod2.pred)))
# examine both mod confusion matrices
mod2.cfmatCV; mod2.cfmat

# calculate accuracies with std.dev = F
mod2.accB <- presence.absence.accuracy(dat2, 
                                       threshold = mod.cut$mod2.pred, 
                                       st.dev = F)
tss <- mod2.accB$sensitivity + mod2.accB$specificity - 1 # code TSS metric
mod2.accB <- cbind(mod2.accB[1:7], tss) # bind all metrics
mod2.accB[c(1, 4:5, 7:8)] # examine accuracies
auc.roc.plot(dat2, color = T)

# spatial probability prediction
spp.predsR <- as(bio_pres_ext, "Raster")
modFprob.LR.1 <- predict(spp.predsR, mod1.LR, filename = "mod2.LRprob.tif", 
                         type = "response", fun = predict, 
                         index = 2, overwrite = T)

# classified prediction
modFprobclas.R=reclassify(modFprob.LR.1,c(0,mod.cut[[2]],0,mod.cut[[2]],1,1))
modFprobclas.bin.R <- reclassify(modFprob.LR.1,c(0,.2,1,
                                                 .2,.4,2,
                                                 .4,.6,3,
                                                 .6,.8,4,
                                                 0.8,1,5))

plot(modFprob.LR.1, main = "Probability Model")
plot(modFprobclas.R, main = "Threshold Classification")
plot(modFprobclas.bin.R, main = "Binned Probability")


# spatial probability prediction - past
lgm.predsR <- as(bio_past_ext, "Raster")
mod.lgmpred.LR.1 <- predict(bio_past_ext, mod1.LR, filename = "mod2.lgmpred.tif", 
                         type = "response", fun = predict, 
                         index = 2, overwrite = T)

# classified prediction
modlgmpredclas.R=reclassify(mod.lgmpred.LR.1,c(0,mod.cut[[2]],0,mod.cut[[2]],1,1))
modlgmpredclas.bin.R <- reclassify(mod.lgmpred.LR.1,c(0,.2,1,
                                                 .2,.4,2,
                                                 .4,.6,3,
                                                 .6,.8,4,
                                                 0.8,1,5))

plot(mod.lgmpred.LR.1, main = "Probability Model")
plot(modlgmpredclas.R, main = "Threshold Classification")
plot(modlgmpredclas.bin.R, main = "Binned Probability")









# now combine into one dataframe
bg.points <- as.data.frame(bg.points)
colnames(bg.points) <- c('Long', 'Lat')
colnames(iv.occ) <- c('Long', 'Lat')
iv.occ$presence <- 1
bg.points$presence <- 0
all_points <- rbind(iv.occ, bg.points)
head(all_points)

### split into testing and training datasets for cross-validation

# we'll arbitrarily set the "testing group" to 1
testing.group <- 1

# create random groups (traditional)
random.group <- kfold(x=all_points, k=10)
table(random.group) # number of observations per group
random.train <- all_points[random.group != testing.group,]
random.test <- all_points[random.group == testing.group,]
parg <- random.train[random.train$presence==1,c('Long','Lat')]
aarg <- random.train[random.train$presence==0,c('Long','Lat')]
############# Fit some distribution models!
# Here we're using Maxent
# fit the model
mod.me <- maxent(x=bioclim_ext, 
                 p=as.data.frame(random.train[random.train$presence==1,c('Long','Lat')]),
                 a=as.data.frame(random.train[random.train$presence==0,c('Long','Lat')]))
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




install.packages("dplyr", version="1.0.10") 

install.packages("gtools", version="3.9.3")

install.packages("plotfunctions", version="1.4")

install.packages("rgbif", version="3.7.3")

install.packages("terra", version="1.6-41")

install.packages("devtools")
library("devtools")
library("dplyr")
library("gtools")
library("plotfunctions")
library("rgbif")
library("terra")
library("megaSDM")
devtools::install_github("brshipley/megaSDM", build_vignettes = TRUE)
setwd("~/yaupon/sdm")

install.packages("dismo")
library("dismo")
install.packages("envirem")
library(envirem)

bradypus <- read_csv(paste0(system.file(package="dismo"), "/ex/bradypus.csv"))
files <- list.files(path=paste(system.file(package="dismo"), '/ex',
                               sep=''), pattern='grd', full.names=TRUE )
# we use the first file to create a RasterLayer
mask <- raster(files[1])


my_extent <- ext(-105, -75, 24, 42) 

# 2. Create the blank template raster object
# Set the resolution (e.g., 0.5 degrees per pixel) and coordinate system (WGS84)
bioclim_ext <- crop(bioclim, my_extent)
terra::plot(my_raster)

chelsa_bio1 <- rast("CHELSA_TraCE21k_bio01_-001_V.1.0.tif")
chelsa_bio1_0020 <- rast("CHELSA_TraCE21k_bio01_0020_V.1.0.tif")
chelsa_bio1_200 <- rast("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1-5/CHELSA_TraCE21k_bio01_-200_V.1.0.tif")
terra::plot(chelsa_bio1)
my_raster <- crop(bioclim[1], my_extent)
mask <- raster(my_raster)
set.seed(1963)
bg <- randomPoints(mask, 500)
par(mfrow=c(1,2))
plot(!is.na(mask), legend=FALSE)
points(bg, cex=0.5)

bio_values_1 = extract(bioclim_ext, bg)
bio_values_df = cbind.data.frame(coordinates(bg), bio_values_1)
head(bio_values_df)

IVV_occurences <- read_csv("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/sdm/occurrences/Ilex_vomitoria.csv")



IVV_data = bioclim %>%
  terra::extract(dplyr::select(IVV_occurences, x, y)) %>%
  bind_cols(IVV_occurences)

logistic_regr_model <- glm(present ~ tmin + precip,
                           family = binomial(link = "logit"),
                           data = hooded_warb_data)
summary(logistic_regr_model)




input_TA <- "~/yaupon/sdm/rasters/wc2-5/*.bil"
RCP45 <- "~/yaupon/sdm/rasters/wc2-5/*.bil"


# If you have your own data, replace the system.file command with
# a pathway to the training area files.

envoutput <- "TestRun"

# Here we define the extent of the training and study regions in c(xmin, xmax, ymin, ymax) form. 
TSEnv <- TrainStudyEnv(input_TA = input_TA,
                       output = envoutput,
                       clipTrain = c(-91.5, -75, 25.5, 36),
                       clipStudy = c(-91.5, -75, 25.5, 36))

Env2050_4.5 <- list.files(system.file("extdata", "predictenv/RCP4.5/2050", package = "megaSDM"), 
                          pattern = ".grd$", 
                          full.names = TRUE)
Env2070_4.5 <- list.files(system.file("extdata", "predictenv/RCP4.5/2070", package = "megaSDM"), 
                          pattern = ".grd$", 
                          full.names = TRUE)
Env4.5 <- list(Env2050_4.5, Env2070_4.5)

# The "time_periods" argument must contain the current time (the time of the training 
# and study rasters) first and then the time periods for the forecast/hindcast.

PredictEnv(studylayers = TSEnv$study,
           futurelayers = Env4.5,
           time_periods = c(2010, 2050, 2070),
           output = envoutput,
           scenario_name = "RCP4.5")

#Repeat with a different climate scenario (RCP8.5):
Env2050_8.5 <- list.files(system.file("extdata", "predictenv/RCP8.5/2050", package = "megaSDM"),
                          pattern = ".grd$", 
                          full.names = TRUE)
Env2070_8.5 <- list.files(system.file("extdata", "predictenv/RCP8.5/2070", package = "megaSDM"), 
                          pattern = ".grd$", 
                          full.names = TRUE)

Env8.5 <- list(Env2050_8.5, Env2070_8.5)
PredictEnv(studylayers = TSEnv$study,
           futurelayers = Env8.5,
           time_periods = c(2010, 2050, 2070),
           output = envoutput,
           scenario_name = "RCP8.5")

# This function only takes occurrences from the described trainingarea extent.
# The defined extent should be the same (or similar to) as the extent of the training area.
# Given in latitude/longitude coordinates:
extent_occ <- c(-91.5, -75, 25.5, 36)

# A list of southeastern mammals for this example
spplist <- "Ilex vomitoria"


# Define the file folder where the occurrences will be written, within the working directory
# (if this folder doesn't already exist, megaSDM will make it)
occ_output <- "occurrences"

Occurrences <- OccurrenceCollection(spplist = spplist,
                                    output = occ_output,
                                    trainingarea = extent_occ)


# NOTE: when running this using R Markdown, you may get "incomplete final line..." 
#    warnings. However, they do not appear to affect the total number or identity
#    of the occurrence points and when the code is run off of the console, the
#    warnings do not appear.

# Because one species was renamed, rename species list to reflect taxonomy changes
spplist <- Occurrences$Scientific.Name


# First, get the list of the occurrence files
occlist <- list.files(occ_output, pattern = ".csv", full.names = TRUE)

OccurrenceManagement(occlist = occlist,
                     output = occ_output,
                     envextract = TRUE,
                     envsample = TRUE,
                     nbins = 25,
                     envdata = TSEnv$training)

# Get the list of occurrence files again, even if they were written out 
# in the same folder as before. This ensures that the occurrence files 
# are properly formatted.  
occlist <- list.files(occ_output, pattern = ".csv", full.names = TRUE)

# The location to print out the background buffers (.shp) (will be created if it doesn't exist)
buff_output <- "TestRun/buffers"

# Generates buffers for each species.
BackgroundBuffers(occlist = occlist,
                  envdata = TSEnv$training,
                  buff_output,
                  ncores = 2)

# Set the parameters for the background point generation 
# (how many points, and how spatially-constrained)

# How many background points should be generated per species?
nbg <- 1000

# What proportion of the background points should be sampled from within the buffers?
spatial_weights <- 0.5 

# Should the background points be environmentally subsampled (Varela) or 
# randomly distributed (random)?
sampleMethod <- "Varela" 

# Because we want a partial spatial constraint (50% of points within the buffer), we must make a 
# list of the buffer files to use in the creation of the background points. In the example, 
# these files are created from the BackgroundBuffers function, but they can also be generated
# outside of megaSDM and brought in here.

bufflist <- list.files(buff_output, pattern = ".shp$", full.names = TRUE)

# Define the location where the background points will be printed out to (.shp) 
# (This directory will be created if it doesn't already exist)
bg_output <- "TestRun/backgrounds"

BackgroundPoints(spplist = spplist,
                 envdata = TSEnv$training,
                 output = bg_output,
                 nbg = nbg,
                 spatial_weights = spatial_weights,
                 buffers = bufflist,
                 method = sampleMethod,
                 ncores = 2)

# Define a list of the environmental variables to keep for each species
# In this example, we simply want all of the species to have the same environmental variables.
envvar <- rep("Bio1,Bio12,Bio14,Bio6,Bio9", length = length(occlist))

# Define a list of the background point files 
# (either created in the BackgroundPoints function or generated separately)
bglist <- list.files(bg_output, pattern = ".csv", full.names = TRUE)

# In this example, megaSDM overwrites the occurrence and background points,
# but they could be placed in a different folder if requested.
VariableEnv(occlist = occlist,
            bglist = bglist,
            env_vars = envvar,
            occ_output = occ_output,
            bg_output = bg_output)

# First, define a list of all background and occurrence point files
occlist <- list.files(occ_output, pattern = ".csv", full.names = TRUE)
bglist <- list.files(bg_output, pattern = ".csv", full.names = TRUE)

# Define where the results of the MaxEnt model runs will be printed out to (as .lambdas files)
model_output <- "TestRun/models"

# "nrep" is set to 4, meaning that the MaxEnt algorithm will run 4 times with different
# subsets of occurrence points for a better representation of the habitat suitability.
MaxEntModel(occlist = occlist,
            bglist = bglist,
            model_output = model_output,
            ncores = 2,
            nrep = 4,
            alloutputs = FALSE)

# First, create a list of the time periods and climate scenarios used in the analysis
# (starting with the year the model is trained on)
time_periods <- c(2010,2050,2070)
scenarios <- c("RCP4.5", "RCP8.5")

# Define the directory where the current study area rasters are located 
#    (generated from the TrainStudyEnv function or brought in from a separate location)
study_dir <- "TestRun/studyarea"

# Define the directories where the future study area rasters are location 
#    (generated from the PredictEnv function or brought in from a separate location)


# Define a list of directories for the projected climate layers, 
# separated into the different climate scenarios and years:
#    list(c(Scenario1Year1, Scenario1Year2), 
#         c(Scenario2Year1, Scenario2Year2))

predictdir <- list(c("TestRun/RCP4.5/2050",
                     "TestRun/RCP4.5/2070"),
                   c("TestRun/RCP8.5/2050",
                     "TestRun/RCP8.5/2070"))

# Define Where the results will be printed out.
# For this example, We'll define a new folder within
# the working directory that is specifically for the model
# projecions and analysis.

result_dir <- "Results"

# Other options are also available (check the documentation page)
MaxEntProj(input = model_output,
           time_periods = time_periods,
           scenarios = scenarios,
           study_dir = study_dir,
           predict_dirs = predictdir,
           output = result_dir,
           aucval = 0.6,
           ncores = 2)

# The time maps will be written out to the directory supplied in "result_dir"
result_dir <- "Results"

createTimeMaps(result_dir = result_dir,
               time_periods = time_periods,
               scenarios = scenarios,
               dispersal = FALSE,
               ncores = 2)
