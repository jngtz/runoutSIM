source("./R/pcm.R")
source("./R/random_walk.R")
source("./R/simulation_to_raster.R")
source("./R/runout_connectivity.R")
source("./R/interactive_plot.R")

# Create three examples
# 1.) using points as source points
# 2.) using top elev %cells as source points
#     ^ can demo determining the simulated connectivity of a runout polygon
# 3.) using source points from a grid

# On git, create folder called DEV, which would have all the functions and 
# the demo script separate... allowing a quick start to manipulate the 
# package without having to depend on building a package...

# Load libraries to handle spatial data ########################################
library(terra)
library(sf)
library(mapview)
library(runoutSim)

# Load data ####################################################################

# Load digital elevation model (DEM)
dem <- rast("Data/elev_nosinks.tif") # use sink filled DEM to remove pits and flats 
# e.g. DMMF::SinkFill(raster::raster()) (our random walk is not an infilling algorithm)

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

# We create a connectivity feature to allow for quicker processing 
# it is essentially an index of all the cells covered by the feature max in 
# the given DEM.

feature_mask <- makeConnFeature(river, dem)

# Run PCM-Random Walk for Single Source Cell ###################################

# Run for a single source point
sim_paths = runoutSim(dem = dem, st_coordinates(source_point), mu = 0.08, md = 140, 
                      slp_thresh = 35, exp_div = 3, per_fct = 1.95, walks = 1000,
                      source_connect = TRUE, connect_feature = feature_mask)

# Convert paths to raster with cell transition frequencies
paths_raster <- walksToRaster(sim_paths, dem)

# Plot results
paths_raster <- crop(paths_raster, ext(runout_polygon)+1000)
plot(paths_raster, legend = T)
plot(st_geometry(river), add = T, border = "#5b86b9")
plot(st_geometry(runout_polygon), add = T)
plot(source_point, add = T)


# Run PCM-Random Walks for Multiple Cells ######################################

# Get coordinates of source points
source_l <- list()
for(i in 1:nrow(source_points)){
  source_l[[i]] <- st_coordinates(source_points[i,])
}

# ^ Need to make a function to clean this up - creating a source object.

# Use lapply to run for multiple source cells
rw_l <- lapply(source_l, function(x) {
  runoutSim(dem = dem, xy = x, mu = 0.08, md = 140, 
            slp_thresh = 50, exp_div = 3.0, per_fct = 1.95, walks = 1000,
            source_connect = TRUE, connect_feature = feature_mask)})

trav_freq <- walksToRaster(rw_l, dem)
trav_prob <- rasterCdf(trav_freq)

# Run Parallel PCM-Random Walks for Multiple Cells #############################

library(parallel)
# Define number of cores to use
n_cores <- detectCores() -2

packed_dem <- wrap(dem)

# Create parallel loop
cl <- makeCluster(n_cores, type = "PSOCK") # Open clusters

# Load objects and 'custom' functions to each cluster
# clusterExport(cl, varlist = c("runoutSim", "euclideanDistance", "adjCells",
#                              "adjRowCol", "pcm", "packed_dem","sourceConnect",
#                              "feature_mask"))

clusterExport(cl, varlist = c("packed_dem", "feature_mask"))

# Load required packages to each cluster
clusterEvalQ(cl, {
  library(terra)
  library(runoutSim)
})


multi_sim_paths <- parLapply(cl, source_l, function(x) {

  runoutSim(dem = unwrap(packed_dem), xy = x, mu = 0.08, md = 140, 
        slp_thresh = 35, exp_div = 3, per_fct = 1.95, walks = 1000,
        source_connect = TRUE, connect_feature = feature_mask)
})

stopCluster(cl) 

trav_freq <- walksToRaster(rw_l, dem)
trav_prob <- rasterCdf(trav_freq)

# Source connectivity ##########################################################

conn_prob <- connToRaster(rw_l, dem)

# Visualize results ############################################################

trav_vel <- velocityToRaster(rw_l, dem)

leafplot(runout_polygons) %>%
  leafplot(trav_freq) %>%
  leafplot(trav_prob) %>%
  leafplot(conn_prob, palette = 'magma') %>%
  leafplot(trav_vel, palette = 'plasma') %>%
  leafplot(source_points, color = "red") %>%
  leafplot(river, color = "#99d2ff")


