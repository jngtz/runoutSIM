# For runout path simulation and optimization
library(runoutSim)
library(runoptGPP)
library(terra)
library(raster)
library(sf)


# Load data ####################################################################

# Load digital elevation model (DEM)
dem <- rast("Dev/Data/elev_fillsinks_WangLiu.tif")

# Compute hillshade for visualization 
slope <- terrain(dem, "slope", unit="radians")
aspect <- terrain(dem, "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270)


# Load debris flow runout source points and polygons
source_points <- st_read("Dev/Data/debris_flow_source_points.shp")
source_points$Id <- NULL
runout_polygons <- st_read("Dev/Data/debris_flow_runout_polygons.shp")



# Using runoptGPP, compute the runout polygon geometry
rungeom <- runoptGPP::runoutGeom(as_Spatial(runout_polygons), raster(dem))
rungeom$fid <- rungeom$id <- NULL
runout_polygons <- cbind(runout_polygons, round(rungeom,2))

# Join runout geometry to corresponding source points
source_points <- st_join(source_points, runout_polygons, join = st_intersects,
                         left = TRUE,
                         largest = FALSE
)

# Load river channel and stream network, and basin boundary
river_channel <- st_read("Dev/Data/river_channel.shp") # may need to buffer to represent high flow conditions...
stream_network <- st_read("Dev/Data/river_rio_olivares.shp")
bnd_catchment <- st_read("Dev/Data/basin_rio_olivares.shp")

# Buffer stream channels
buffer_stream <- st_buffer(stream_network, dist = 30)

# Combine river and stream network into one polygon
drainage_network <- st_sf(st_union(st_union(st_geometry(river_channel), st_geometry(buffer_stream), is_coverage = TRUE)))

# Create an interactive leaflet map to explore data
leafmap(bnd_catchment, color = '#f7f9f9', fill_color = '#FF000000', weight = 4) %>%
    leafmap(runout_polygons) %>%
    leafmap(source_points, color = '#e74c3c') %>%
    leafmap(hill, palette = grey(0:100/100), opacity = 1, add_legend = FALSE) %>%
    leafmap(dem, palette = viridis::mako(10), opacity = 0.6) %>%
    leafmap(drainage_network, fill_color = "lightblue") 


# Simulate runout ##############################################################

# Get coordinates of source points to create a source list object
source_l <- makeSourceList(source_points)

# Perform random walk simulation for each source point
sim_runs <- list()

for(i in 1:length(source_l)){
  
  print(paste(i, "of", length(source_l), ":", Sys.time()))
  
  sim_runs[[i]] <- runoutSim(dem = dem, xy = source_l[[i]], 
                              mu = 0.08, 
                              md = 40, 
                              slp_thresh = 40, 
                              exp_div = 3, 
                              per_fct = 1.9, 
                              walks = 1000)
}


# Coerce results to a raster
trav_freq <- walksToRaster(sim_runs, method = "freq", dem)
vel_ms <- velocityToRaster(sim_runs, dem)
trav_prob <- runoutSim::rasterCdf(walksToRaster(sim_runs, method = "freq", dem))


# Map / visualize the results ##################################################

leafmap(bnd_catchment, color = 'black', fill_color = '#FF000000', weight = 4,
        opacity = 1, label = "Basin boundary") %>%
    # Debris flow observations
  leafmap(runout_polygons, opacity = 0.3, label = "Runout polygons") %>%
  leafmap(source_points, color = '#e74c3c', label = "Source points") %>%
  # DEM data
  leafmap(hill, palette = grey(0:100/100), opacity = 1, add_legend = FALSE,
          label = "Hillshade") %>%
  leafmap(dem, palette = viridis::mako(10), opacity = 0.6, label = "DEM") %>%
  # Modelling results
  leafmap(trav_prob,  palette = viridis::viridis(10, direction = -1), 
          label = "Trav. Prob.") %>%
  leafmap(trav_freq,  palette = viridis::viridis(10, direction = -1),
          label = "Trav. Freq.") %>%
  leafmap(vel_ms,  palette = viridis::viridis(10, direction = -1),
          label = "Vel. ms^-1")
