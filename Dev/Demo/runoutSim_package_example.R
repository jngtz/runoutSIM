# Create three examples
# 1.) using points as source points
# 2.) using top elev %cells as source points
#     ^ can demo determining the simulated connectivity of a runout polygon
# 3.) using source points from a grid

# On git, create folder called DEV, which would have all the functions and 
# the demo script separate... allowing a quick start to manipulate the 
# package without having to depend on building a package...

# Load libraries to handle spatial data ########################################
library(runoutSim)
library(terra)
library(sf)

# Load data ####################################################################

# Load digital elevation model (DEM)
dem <- rast("Dev/Data/elev_fillsinks_WangLiu.tif")# use sink filled DEM to remove pits and flats 
dem <- rast("Dev/Data/elev.tif")
# e.g. DMMF::SinkFill(raster::raster()) (our random walk is not an infilling algorithm)

# Hillshade for visualization 
slope <- terrain(dem, "slope", unit="radians")
aspect <- terrain(dem, "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270)

# Load runout source points and polygons
source_points <- st_read("Dev/Data/debris_flow_source_points.shp")
source_points$run_id <- 1:nrow(source_points)

runout_polygons <- st_make_valid(st_read("Dev/Data/debris_flow_runout_polygons.shp"))
runout_polygons$run_id <- 1:nrow(runout_polygons)
# ^ Need to clean this up so make valid not needed

# Select a single debris flow and source point for the example
runout_polygon <- runout_polygons[50,]

# Get corresponding source point
source_point  <- st_filter(st_as_sf(source_points), st_as_sf(runout_polygon))

# Load river data used for connecivity analysis

river_channel <- st_read("Dev/Data/river_channel.shp")
stream_channels <- st_read("Dev/Data/river_rio_olivares.shp")

bnd_catchment <- st_read("Dev/Data/basin_rio_olivares.shp")
# buffer stream channels
buffer_stream <- st_buffer(stream_channels, dist = 30)

# combine to one feature
drainage_network <- st_sf(st_union(st_union(st_geometry(river_channel), st_geometry(buffer_stream), is_coverage = TRUE)))


# rasterize drainage_network to view cells that will be considered for connectivity
r_dn = rasterize(drainage_network, dem)

# View data
leafmap(bnd_catchment, 
           color = '#f7f9f9', 
           fill_color = '#FF000000', 
           weight = 4) %>%
  
  leafmap(runout_polygons) %>%
  
  leafmap(source_points, 
           color = '#e74c3c') %>%
  
  
  leafmap(hill, palette = grey(0:100/100), opacity = 1, add_legend = FALSE) %>%
  
  leafmap(dem, palette = viridis::mako(10), opacity = 0.6) %>%
  
  leafmap(r_dn, 
           label = "Drainage Network", 
           palette = list(classes = 1, colors = "#99d2ff", labels = ""), 
           opacity = 0.8)
  

# Data pre-processing ##########################################################

# We create a connectivity feature to allow for quicker processing 
# it is essentially an index of all the cells covered by the feature max in 
# the given DEM.

feature_mask <- makeConnFeature(drainage_network, dem)

# Run PCM-Random Walk for Single Source Cell ###################################

# Run for a single source point
sim_paths = runoutSim(dem = dem, st_coordinates(source_point), mu = 0.49, md = 40, 
                      slp_thresh = 40, exp_div = 3, per_fct = 1.9, walks = 1000,
                      source_connect = TRUE, connect_feature = feature_mask)

# Convert paths to raster with cell transition frequencies
paths_raster <- walksToRaster(sim_paths, dem)

# Plot results
paths_raster <- crop(paths_raster, ext(runout_polygon)+1000)
plot(paths_raster, legend = T)
#plot(st_geometry(drainage_network), add = T, border = "#5b86b9")
plot(st_geometry(runout_polygon), add = T)
plot(st_geometry(source_point), add = T)


# Run PCM-Random Walks for Multiple Cells ######################################

# Get coordinates of source points
source_l <- list()
for(i in 1:nrow(source_points)){
  source_l[[i]] <- st_coordinates(source_points[i,])
}

# ^ Need to make a function to clean this up - creating a source object.
for(i in 1:length(source_l)){
  print(paste("source", i))
  runoutSim(dem = dem, xy = source_l[[i]], mu = 0.08, md = 140, 
            slp_thresh = 50, exp_div = 3.0, per_fct = 1.95, walks = 1000,
            source_connect = TRUE, connect_feature = feature_mask)
}


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

trav_freq <- walksToRaster(multi_sim_paths, dem)
trav_prob <- rasterCdf(trav_freq)

# Source connectivity ##########################################################

conn_prob <- connToRaster(multi_sim_paths, dem)

# Visualize results ############################################################

trav_vel <- velocityToRaster(multi_sim_paths, dem)

leafmap(bnd_catchment, 
         color = '#2e4053', 
         fill_color = '#FF000000', 
         weight = 4) %>%
  leafmap(drainage_network, color = "#99d2ff") %>%
  leafmap(runout_polygons) %>%
  leafmap(trav_freq) %>%
  leafmap(trav_prob) %>%
  leafmap(trav_vel, palette = 'plasma') %>%
  leafmap(conn_prob, palette = 'magma') %>%
  leafmap(source_points, color = "red")


# add export as html...

