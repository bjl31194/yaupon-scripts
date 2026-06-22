
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
