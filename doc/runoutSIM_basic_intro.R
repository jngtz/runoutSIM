## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.show='hold'----------------------------------------------------------
# load packages
library(runoutSim)
library(terra)
library(sf)

# Load digital elevation model (DEM)
dem <- rast("C:/GitProjects/runoutSIM/Data/elev_fillsinks_WangLiu.tif")

# Compute hillshade for visualization 
slope <- terrain(dem, "slope", unit="radians")
aspect <- terrain(dem, "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270)

# Load debris flow runout source points and polygons
source_points <- st_read("C:/GitProjects/runoutSIM/Data/debris_flow_source_points.shp")
runout_polygons <- st_read("C:/GitProjects/runoutSIM/Data/debris_flow_runout_polygons.shp")

# Load basin boundary
bnd_catchment <- st_read("C:/GitProjects/runoutSIM/Data/basin_rio_olivares.shp")


## ----fig.show='hold', out.width="100%", out.height="500px"--------------------
map <- leafmap(bnd_catchment, color = '#f7f9f9', fill_color = '#FF000000', weight = 4) %>%
    leafmap(runout_polygons) %>%
    leafmap(source_points, color = '#e74c3c') %>%
    leafmap(hill, palette = grey(0:100/100), opacity = 1, add_legend = FALSE, 
            add_image_query = FALSE) %>%
    leafmap(dem, palette = viridis::mako(10), opacity = 0.6, add_image_query = FALSE)
map

## -----------------------------------------------------------------------------
# Select a single debris flow and source point for the example
runout_polygon <- runout_polygons[31,]

# Get corresponding source point
source_point  <- st_filter(st_as_sf(source_points), st_as_sf(runout_polygon))

# Get coordinates of source point
source_crds <- st_coordinates(source_point)
print(source_crds)

## ----fig.show='hold'----------------------------------------------------------
sim_paths = runoutSim(dem = dem, 
                      xy = source_crds, 
                      mu = 0.08, 
                      md = 140, 
                      slp_thresh = 35, 
                      exp_div = 2.5, 
                      per_fct = 1.95, 
                      walks = 1000)

# plot structure of sim_paths
str(sim_paths)

## ----fig.height=5, fig.width = 7, fig.align='center'--------------------------
## Convert paths to raster with cell transition frequencies
paths_raster <- walksToRaster(sim_paths, dem, method = "freq")

# Plot results
paths_raster <- crop(paths_raster, ext(runout_polygon)+500)
plot(crop(hill, ext(runout_polygon)+500), col=grey(0:100/100), legend=FALSE,
     mar=c(2,2,1,4), main = "Runout traverse frequency")
plot(paths_raster, legend = T, alpha = 0.5, add = TRUE)
plot(st_geometry(runout_polygon), add = TRUE)
plot(source_point, add = TRUE)

## ----fig.height=5, fig.width = 7, fig.align='center'--------------------------
## Convert paths to raster with max. velocities
vel_raster <- velocityToRaster(sim_paths, dem)

# Plot results
vel_raster <- crop(vel_raster, ext(runout_polygon)+500)
plot(crop(hill, ext(runout_polygon)+500), col=grey(0:100/100), legend=FALSE,
     mar=c(2,2,1,4), main = "Runout velocity")
plot(vel_raster, col = map.pal('plasma'), legend = TRUE, alpha = 0.5, add = TRUE)
plot(st_geometry(runout_polygon), add = TRUE)
plot(source_point, add = TRUE)

## -----------------------------------------------------------------------------
# Load river channel polygon (vector)
river_channel <- st_read("C:/GitProjects/runoutSIM/Data/river_channel.shp")

# Create a connectivity feature for runoutSim that matches the input DEM
feature_mask <- makeConnFeature(river_channel, dem)

## ----fig.height=7, fig.width = 5, fig.align='center'--------------------------
sim_paths = runoutSim(dem = dem, 
                      xy = source_crds, 
                      mu = 0.08, 
                      md = 140, 
                      slp_thresh = 35, 
                      exp_div = 2.5, 
                      per_fct = 1.95, 
                      walks = 1000,
                      source_connect = TRUE,
                      connect_feature = feature_mask)

# Plot structure of sim_paths
str(sim_paths)

# Get connectivity probability
sim_paths$prob_connect

tp_raster <- walksToRaster(sim_paths, dem, method = "cdf_prob")
conn_raster <- connToRaster(sim_paths, dem)

# Crop results
tp_crop <- crop(tp_raster, ext(runout_polygon)+ 500)
conn_crop <- crop(conn_raster, ext(runout_polygon)+ 500)

# Plot results
par(mfrow = c(2,1))
plot(crop(hill, ext(runout_polygon)+500), col=grey(0:100/100), legend=FALSE,
     mar=c(2,2,1,4), main = "Traverse prob.")
plot(st_geometry(river_channel), col = "lightblue", add = TRUE)
plot(tp_crop, legend = TRUE, alpha = 0.5, add = TRUE)
plot(st_geometry(runout_polygon), add = TRUE)
plot(source_point, add = TRUE)


plot(crop(hill, ext(runout_polygon)+500), col=grey(0:100/100), legend=FALSE,
     mar=c(2,2,1,4), main = "Connectivity prob.")
plot(st_geometry(river_channel), col = "lightblue", add = TRUE)
plot(conn_crop, legend = TRUE, alpha = 0.5, add = TRUE)
plot(st_geometry(runout_polygon), add = TRUE)
plot(source_point, add = TRUE)

## ----fig.show='hold', out.width="100%", out.height="500px"--------------------
# Get coordinates of source points to create a source list object
source_l <- makeSourceList(source_points)

# Perform random walk simulation for each source point
sim_runs <- list()

for(i in 1:length(source_l)){
  
  sim_runs[[i]] <- runoutSim(dem = dem, xy = source_l[[i]], 
                              mu = 0.08, 
                              md = 40, 
                              slp_thresh = 40, 
                              exp_div = 3, 
                              per_fct = 1.9, 
                              walks = 1000)
}




## ----fig.show='hold', fig.height=4, fig.width = 7, fig.align='center'---------
# Coerce results to a raster
trav_freq <- walksToRaster(sim_runs, method = "freq", dem)
vel_ms <- velocityToRaster(sim_runs, dem, method = "max")
trav_prob <- walksToRaster(sim_runs, method = "max_cdf_prob", dem)


# Plot random walks from mulitple source cells
par(mfrow = c(1,3))
plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4), 
     main = "Traverse frequency")
plot(trav_freq, add = TRUE, alpha = 0.7)

plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4), 
     main = "Traverse probability (max.)")
plot(trav_prob, add = TRUE, alpha = 0.7)

plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4), 
     main = "Velocity (m/s)")
plot(vel_ms, col = map.pal("plasma"), add = TRUE, alpha = 0.7)



## ----fig.show='hold', out.width="100%", out.height="500px"--------------------

# Create interactive leaflet map
sim_map <- 
  # Debris flow observations
  leafmap(runout_polygons, opacity = 0.3, label = "Runout polygons") %>%
  leafmap(source_points, color = '#e74c3c', label = "Source points") %>%
  # Terrain data
  leafmap(hill, palette = grey(0:100/100), opacity = 1, add_legend = FALSE, 
            add_image_query = FALSE, label = "Hillshade") %>%
  leafmap(dem, palette = viridis::mako(10), opacity = 0.6, 
          add_image_query = FALSE, label = "DEM") %>%
  # Modelling results
  leafmap(trav_prob,  palette = viridis::viridis(10, direction = -1), 
          label = "Traverse prob", add_image_query = FALSE) %>%
  leafmap(vel_ms,  palette = viridis::plasma(10, direction = -1),
          label = "Velocity", add_image_query = FALSE)

# Start with mapped zoomed in
sim_map <- leaflet::setView(sim_map, lng = -70.13694, lat = -33.3703, zoom = 13)


# Hide layers
sim_map <- leaflet::hideGroup(sim_map, c("Velocity", "Source points", "DEM", "Hillshade"))

sim_map


