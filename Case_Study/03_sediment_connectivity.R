
# Load libraries to handle spatial data ########################################
library(runoutSim)
library(terra)
library(sf)

# Load data ####################################################################

# Load digital elevation model (DEM)
dem <- rast("Dev/Data/elev_nosinks.tif") # use sink filled DEM to remove pits and flats 
# e.g. DMMF::SinkFill(raster::raster()) (our random walk is not an infilling algorithm)
source_areas <- rast("Data/auto_classified_source_areas.tif")
# Hillshade for visualization 
slope <- terrain(dem, "slope", unit="radians")
aspect <- terrain(dem, "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270)

# Load runout source points and polygons
source_points <- st_read("Dev/Data/debris_flow_source_points.shp")
runout_polygons <- st_make_valid(st_read("Dev/Data/debris_flow_runout_polygons.shp"))
# ^ Need to clean this up so make valid not needed

# Select a single debris flow and source point for the example
runout_polygon <- runout_polygons[10,]

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
          opacity = 0.8) %>%
  
  leafmap(source_areas, palette = list(classes = 1, colors = "#e5207a", labels = "Source area"))


# We create a connectivity feature to allow for quicker processing 
# it is essentially an index of all the cells covered by the feature max in 
# the given DEM.

feature_mask <- makeConnFeature(drainage_network, dem)

source_areas<- rast("Data/auto_classified_source_areas.tif") # from "02_source_area_prediction.R"
source_areas <- crop(source_areas, dem)
sum(values(source_areas), na.rm = TRUE)

# Find cells where the value is 1
source_cells <- which(values(source_areas) == 1)

# Extract the coordinates of these cells
source_xy <- xyFromCell(source_areas, source_cells)

# Let's do a random sub-sample for testing
#test_sample_index <- sample(1:nrow(source_xy), size = 1000)
#source_xy <- source_xy[1:100,]

# Create a list of matricies for input to pcmRW
source_l <- makeSourceList(source_xy)



library(parallel)

setwd("C:\\sda\\Workspace\\sedconnect")

(load("Global_MBO_RunoutSim.Rd"))

packed_dem <- wrap(dem)

n_cores <- detectCores() -2
cl <- makeCluster(n_cores) # Open clusters
clusterExport(cl, varlist = c("global_run", "packed_dem", "feature_mask"))

# Load required packages to each cluster
clusterEvalQ(cl, {
  library(terra);
  library(runoutSim)
})

print(start_time <- Sys.time())

multi_sim_paths <- parLapply(cl, source_l, function(x) {
  
  local_dem <- terra::unwrap(packed_dem)
  
  runoutSim(dem = local_dem, xy = x, 
            mu = global_run$x$mu, 
            md = global_run$x$md, 
            slp_thresh = global_run$x$slp, 
            exp_div = global_run$x$ex, 
            per_fct = global_run$x$per, 
            walks = 1000,
            connect_feature = feature_mask,
            source_connect = TRUE)
})
print(end_time <- Sys.time())
print(run_time <- end_time - start_time)



save(multi_sim_paths, file = "runoutSim_wConnectFeature.Rd")
stopCluster(cl) # Close clusters


for(i in 1:10){
  print(i)
  local_dem <- terra::unwrap(packed_dem)
  
  runoutSim(dem = local_dem, xy = source_l[[i]], 
            mu = 0.14, 
            md = 40, 
            slp_thresh = 30, 
            exp_div = 2, 
            per_fct = 1.5, 
            walks = 1000,
            connect_feature = feature_mask,
            source_connect = FALSE)
}