# runoutSIM case study: Source area prediction #################################

# This script is 2 of 3 scripts presenting a case study (demo) for the runoutSIM
# package for regional mass movement runout simulation.

# It performs the spatial classification of potential source areas (cells) using
# statistical learning.

# By Jason Goetz, PhD (jgoetz@wlu.ca)
# Department of Geography and Environmental Studies,
# Wilfrid Laurier University, Canada

# August 7, 2025

## Steps: ####
# 1. Load required packages
# 2. Load required spatial data
# 3. Generate sample data
# 4. Build spatial prediction model
# 5. Model validation
# 6. Apply spatial prediction
# 7. Classify prediction map

# Load required packages #######################################################

library(terra)
library(sp) 
library(sf)
library(runoutSim)


# Load data ####################################################################

# Read mapped source points and runout polygon shapefiles
source_points <- st_read("Data/debris_flow_source_points.shp")
runout_polygons <- st_read("Data/debris_flow_runout_polygons.shp")

# Read river channel and catchment boundary shapefile
river_channel <- st_read("Data/river_channel.shp")

bnd_catchment <- st_read("Data/basin_rio_olivares.shp")

## Visualization input spatial data ####
leafmap(source_points, color = "red") %>% leafmap(runout_polygons) %>% 
  leafmap(river_channel, col = "#2e86c1") %>% 
  leafmap(bnd_catchment, col = "white", fill_color = "#FF000000",
          weight = 4)

river_channel <- st_read("Data/river_channel.shp") # may need to buffer to represent high flow conditions...

# Read elevation model raster using raster package
dem <- rast("Data/elev_fillsinks_WangLiu.tif")
dem

# Note: This DEM is from the (free) publicly available ALOS PALSAR Radiometric Terrain Corrected
#       (RTC) high-resolution (12.5 m) digital elevation model.

#       ASF DAAC: ALOS PALSAR Radiometric Terrain Corrected high
#         resolution digital elevation model. Includes Material ©JAXA/METI
#         (2011). DOI: 10.5067/Z97HFCNKR6VA.

#       https://asf.alaska.edu/data-sets/derived-data-sets/alos-palsar-rtc/alos-palsar-radiometric-terrain-correction/


# Read terrain attribute data
layer_files <- c("Data/elev.tif",
                 "Data/cslope.tif",
                 "Data/slope.tif",
                 "Data/logcarea.tif",
                 "Data/plan_curv.tif",
                 "Data/prof_curv.tif",
                 "Data/swi.tif")

layers <- rast(layer_files)
plot(layers)

# Note: These terrain attributes were processed using SAGA-GIS.
#       First, a sink filled DEM was computed using the "Fill Sinks
#       (Wang and Liu, 2006)" tool. Using this sink-filled DEM,
#       the "SAGA Wetness Index" and the "Slope, Aspect, Curvature"
#       tools were used to compute the terrain attributes.

# Create mask of study area
mask_v <- values(layers$elev)
mask_v[!is.na(mask_v)] = 1
mask_area <- setValues(layers$elev, mask_v)

# Mask out runout polygons as well
mask_area <- mask(mask_area,vect(runout_polygons), inverse = TRUE)
mask_area <- mask(mask_area, vect(river_channel), inverse = TRUE)

# Create 
runout_area <- setValues(layers$elev, mask_v)
runout_area <- mask(runout_area, vect(runout_polygons))
runout_area <- mask(runout_area, vect(river_channel), inverse = TRUE)
plot(runout_area)


# Generate sample data #########################################################

# To provide data for our spatially predictive models we need to create
# sample of the response variable: slide and non-slides.

# Set the seed of R's random number generator to reproduce the random samples
# if the script is run at a later time
set.seed(1234)

## Generate random samples ####
# First, we randomly sample cell locations within the runout polygons
smp.slides <- spatSample(runout_area, size = 500, xy = TRUE, as.df = TRUE,
                         method = "random", na.rm = TRUE)

# Remove raster value attribute
smp.slides[,3] = NULL

plot(runout_area)
points(smp.slides)

# Randomly sample cell locations of non-slide cells
smp.noslides <- spatSample(mask_area, size = 500, xy = TRUE, as.df = TRUE,
                           method = "random", na.rm = TRUE)
smp.noslides[,3] = NULL
points(smp.noslides, pch = 20)

## Extract predictor variable values ####

# Extract predictor variable values from grids for the slide samples
df.slides <- extract(layers, smp.slides, sp = TRUE)
df.slides <- as.data.frame(df.slides)
# add xy values
df.slides$x <- smp.slides$x
df.slides$y <- smp.slides$y
# add variable indicating that these are slides
df.slides$slide <- "TRUE"

# do the same for the non-slide sample
df.noslides <- as.data.frame(extract(layers, smp.noslides, sp = TRUE))
df.noslides$x <- smp.noslides$x
df.noslides$y <- smp.noslides$y
df.noslides$slide <- "FALSE"

# combine the slides and no slides data into one dataframe
d <- rbind(df.slides, df.noslides)

# encode the slide variable as a factor
d$slide <- as.factor(d$slide)


# Build spatial prediction model ###############################################

# Load `mgcv` package for generalized additive modelling
library(mgcv)

# make model formula
fml.gam <- slide ~ s(elev, k = 4) + s(slope, k = 4) + s(cslope, k = 4) +
  s(plan_curv, k = 4) + s(prof_curv, k = 4) + s(logcarea, k = 4) +
  s(swi, k = 4)

# fit a generalized additive model (GAM)
model.gam <- gam(fml.gam, family = binomial, data = d)

# visualize smoothing functions in the GAM
par(mfrow=c(2,4))
plot(model.gam)
par(mfrow=c(1,1))

# Model validation ############################################################

# Load library for spatial error estimation
library(sperrorest)

# Perform 5-fold 10-repeated spatial cross-validation
gam.results <- sperrorest(formula = fml.gam, data = d, coords = c("x","y"),
                          model_fun = gam,
                          model_args = list(family=binomial),
                          pred_args = list(type="response"),
                          smp_fun = partition_kmeans,
                          smp_args = list(repetition = 1:10, nfold = 5))

# Get the area under the receiver operating characteristic curve (AUROC)
# values from spatial CV results

gam.auroc.test <- unlist(summary(gam.results$error_rep,level=1)[,"test_auroc"])

# Get average and standard deviation performance value across repititions and
# cross-validation iterations

mean(gam.auroc.test)
sd(gam.auroc.test)


# Apply spatial predictions  ###################################################

# Transform any (raster) layers to match model input values
layers$plan_curv[layers$plan_curv > 0.1] <- 0.1
layers$plan_curv[layers$plan_curv < -0.1] <- -0.1

layers$prof_curv[layers$plan_curv > 0.1] <- 0.1
layers$prof_curv[layers$plan_curv < -0.1] <- -0.1

# Apply prediction to raster data
pred.gam<- terra::predict(layers, model.gam, type = "response", progress = TRUE)

# Mask out areas in river channel bead
maskpred.gam <- mask(pred.gam, vect(river_channel), inverse = TRUE)

# export to a raster format
writeRaster(maskpred.gam, "Data/gam_src_pred_mask.tif", overwrite = TRUE)

maskpred.gam <- rast("Data/gam_src_pred_mask.tif")

## Visualize prediction map #####

plot(maskpred.gam)
leafmap(maskpred.gam, label = "Source area probability")


# Classify prediction map ######################################################

# Apply a typical class membership threshold of 0.5 to classify grid cells as 
# being a source (release) cell

src_class <- maskpred.gam

src_class[src_class >= 0.5] <- 1
src_class[src_class < 0.5] <- 0
src_class[is.na(src_class)] <- 0

plot(src_class)

# The classification contains some issolated (noisy) grid cells. We will apply
# a post-classification filter (majority) to remove them.

# Define a majority filter function (no default in R)
majority <- function(x) {
  x <- x[!is.na(x)]  # Remove NAs
  if (length(x) == 0) return(NA)  # Handle empty case
  return(as.numeric(names(which.max(table(x)))))  # Most frequent value
}

# Apply a majority filter with a 7x7 (rectangular) moving window
x_filter <- focal(src_class, w=7, fun=majority, na.policy="all")
plot(x_filter)

# Mask using boundaries of DEM
mx_filter <- mask(x_filter, dem)
mx_filter[mx_filter == 0] <- NA

## Visualize post-classification

leafmap(mx_filter, palette = list(classes = 1, colors = "#e5207a", 
                                  labels = "Source area"))

## Export source area classification ####
writeRaster(mx_filter, filename = "Data/classified_w7filter_source_areas.tif",
            overwrite = TRUE)
