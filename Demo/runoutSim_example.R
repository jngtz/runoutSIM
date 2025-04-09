source("./R/pcm.R")
source("./R/random_walk.R")
source("./R/combine_walks.R")
source("./R/runout_connectivity.R")
source("./R/interactive_plot.R")

# Load libraries to handle spatial data ########################################
library(terra)
library(sf)
library(mapview)

# Load data ####################################################################

# Load digital elevation model (DEM)
dem <- rast("Data/elev.tif")

# Load runout source points and polygons
source_points <- st_read("Data/debris_flow_source_points.shp")
runout_polygons <- st_make_valid(st_read("Data/debris_flow_runout_polygons.shp"))
# ^ Need to clean this up so make valid not needed

# Select a single debris flow and source point for the example
runout_polygon <- runout_polygons[10,]

# Get corresponding source point
source_point  <- st_filter(st_as_sf(source_points), st_as_sf(runout_polygon))

# Load river data used for connecivity analysis
river <- st_read("Data/river_channel.shp")

# Data pre-processing ##########################################################

# Create a feature for connectivity (object) analysis using the river
# ^ turn this into a function

makeConnFeature <- function(x,y){
  object = terra::rasterize(vect(x), y)
  
  # Create feature mask outside of function (just raster converted to matrix)
  feature_mask <- as.matrix(object, wide=TRUE)
  feature_mask[is.na(feature_mask)] <- 0
  return(feature_mask)
}

feature_mask <- makeConnFeature(river, dem)

# Run PCM-Random Walk for Single Source Cell ##########################################################

# Run for a single source point
sim_paths = runoutSim(dem = dem, st_coordinates(source_point), mu = 0.08, md = 140, 
                      slp_thresh = 35, exp_div = 3, per_fct = 1.95, walks = 1000,
                      source_connect = TRUE, feature_layer = feature_mask)

# Convert paths to raster with cell transition frequencies
paths_raster <- walksToRaster(sim_paths, dem)

# Plot results
paths_raster <- crop(paths_raster, ext(runout_polygon)+1000)
plot(paths_raster, legend = T)
plot(st_geometry(river), add = T, border = "#5b86b9")
plot(st_geometry(runout_polygon), add = T)
plot(source_point, add = T)



