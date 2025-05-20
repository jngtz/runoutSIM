# Load libraries to handle spatial data ########################################
library(runoutSim)
library(terra)
library(sf)

# Load data ####################################################################

# Load digital elevation model (DEM)
dem <- rast("Dev/Data/elev_nosinks.tif") # use sink filled DEM to remove pits and flats 
# e.g. DMMF::SinkFill(raster::raster()) (our random walk is not an infilling algorithm)

# Hillshade for visualization 
slope <- terrain(dem, "slope", unit="radians")
aspect <- terrain(dem, "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270)

# Load runout source points and polygons
source_points <- st_read("Dev/Data/debris_flow_source_points.shp")
runout_polygons <- st_make_valid(st_read("Dev/Data/debris_flow_runout_polygons.shp"))
runout_polygons$sim_id <- 1:nrow(runout_polygons)
# ^ Need to clean this up so make valid not needed

# Select a single debris flow and source point for the example
runout_polygon <- runout_polygons[33,]

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


sim_paths = runoutSim(dem = dem, st_coordinates(source_point), mu = 0.08, md = 140, 
                      slp_thresh = 35, exp_div = 3, per_fct = 1.95, walks = 1000)

# Convert paths to raster with cell transition frequencies
paths_raster <- walksToRaster(sim_paths, dem)

# Plot results
paths_raster <- crop(paths_raster, ext(runout_polygon)+500)
plot(paths_raster, legend = T)
plot(st_geometry(drainage_network), add = T, border = "#5b86b9")
plot(st_geometry(runout_polygon), add = T)
plot(source_point, add = T)

par(mfrow = c(2,1))
# Traverse prob
trav_prob <- walksToRaster(sim_paths, dem) / 1000
trav_prob <- crop(trav_prob, ext(runout_polygon)+300)
hill_crop <- crop(hill, trav_prob)
#plot(hill_crop, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(trav_prob, legend = T, axes = F)
plot(source_point, pch = 19, col = "black", add = T)

trav_ecdf <- rasterCdf(walksToRaster(sim_paths, dem))
trav_ecdf <- crop(trav_ecdf, ext(runout_polygon)+300)
hill_crop <- crop(hill, trav_ecdf)
#plot(hill_crop, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(trav_ecdf, legend = T, axes = F)
plot(source_point, pch = 19, col = "black", add = T)



library(terra)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)

# Convert rasters to data frames
df_freq <- as.data.frame(paths_raster, xy = TRUE, na.rm = TRUE)
names(df_prob)[3] <- "freq"

df_ecdf <- as.data.frame(trav_ecdf, xy = TRUE, na.rm = TRUE)
names(df_ecdf)[3] <- "ecdf"

# Convert source point to data frame
source_df <- as.data.frame(source_point, xy = TRUE)

# Plot traversal probability
m.freq <- ggplot(df_freq, aes(x = x, y = y, fill = freq)) +
  geom_raster() +
  scale_fill_viridis(name = "Traverse\nFrequency", na.value = "transparent") +
  coord_equal() +
  theme_void() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.position = "bottom", legend.direction = "horizontal"
  ) +
  xlab("") + ylab("") +
  coord_fixed() 

# Plot ECDF
m.ecdf <- ggplot(df_ecdf, aes(x = x, y = y, fill = ecdf)) +
    geom_raster() +
    scale_fill_viridis(name = "Traverse\nProbability\n(ECDF)", na.value = "transparent") +
    coord_equal() +
    theme_void() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.position = "bottom", legend.direction = "horizontal"
  ) +
  xlab("") + ylab("") +
  coord_fixed() 

  
m.trav_probs <- m.freq + m.ecdf

ggsave("Case_Study/Figures/traversal_maps.png", plot = m.trav_probs, width = 170, height = 85, units = "mm", dpi = 300)
