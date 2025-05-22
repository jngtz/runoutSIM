# Load libraries to handle spatial data ########################################
library(runoutSim)
library(terra)
library(sf)

# Load data ####################################################################

# Load digital elevation model (DEM)
dem <- rast("Dev/Data/elev_nosinks.tif") # use sink filled DEM to remove pits and flats 

src_pred <- rast("Data/src_pred_mask.tif")
src_area <- rast("Data/auto_classified_source_areas.tif")
# e.g. DMMF::SinkFill(raster::raster()) (our random walk is not an infilling algorithm)

# Hillshade for visualization 
slope <- terrain(dem, "slope", unit="radians")
aspect <- terrain(dem, "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270)

bnd_catchment <- st_read("Dev/Data/basin_rio_olivares.shp")
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

run_vel <- velocityToRaster(sim_paths, dem)
run_vel <- crop(run_vel, ext(runout_polygon)+300)

plot(run_vel, legend = T, axes = F)
plot(source_point, pch = 19, col = "black", add = T)


library(terra)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)

# Convert rasters to data frames
df_freq <- as.data.frame(paths_raster, xy = TRUE, na.rm = TRUE)
names(df_freq)[3] <- "freq"

df_ecdf <- as.data.frame(trav_ecdf, xy = TRUE, na.rm = TRUE)
names(df_ecdf)[3] <- "ecdf"

df_vel <- as.data.frame(run_vel, xy = TRUE, na.rm = TRUE)

# Convert source point to data frame
source_df <- as.data.frame(source_point, xy = TRUE)

# Plot traversal probability
m.freq <- ggplot(df_freq, aes(x = x, y = y, fill = freq)) +
  geom_raster() +
  scale_fill_viridis(name = "Traverse\nFrequency", na.value = "transparent",
                     direction = 1) +
  coord_equal() +
  theme_void() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.position = "bottom", legend.direction = "horizontal"
  ) +
  xlab("") + ylab("") +
  coord_sf() +
  annotation_scale(aes(style = "ticks", location = "tr"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) 

# Plot ECDF
m.ecdf <- ggplot(df_ecdf, aes(x = x, y = y, fill = ecdf)) +
    geom_raster() +
    scale_fill_viridis(name = "Traverse\nProbability\n(ECDF)", 
                       na.value = "transparent",
                     direction = 1) +
    coord_equal() +
    theme_void() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.position = "bottom", legend.direction = "horizontal"
  ) +
  xlab("") + ylab("") +
  coord_sf() + 
  annotation_scale(aes(style = "ticks", location = "tr"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) 

m.vel <- ggplot(df_vel, aes(x = x, y = y, fill = max_velocity_ms)) +
  geom_raster() +
  scale_fill_viridis(name = "Maximum\nVelocity\n(m/s)", 
                     option = "plasma", 
                     direction =  1, 
                     na.value = "transparent") +
  coord_equal() +
  theme_void() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.position = "bottom", legend.direction = "horizontal"
  ) +
  xlab("") + ylab("") +
  coord_sf() + 
  annotation_scale(aes(style = "ticks", location = "tr"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) 

  
m.trav_probs <- m.freq + m.ecdf

m.output_single <- m.freq + m.ecdf + m.vel


ggsave("Case_Study/Figures/traversal_maps.png", plot = m.trav_probs, width = 170, height = 85, units = "mm", dpi = 300)
ggsave("Case_Study/Figures/eg_single_maps.png", plot = m.output_single, width = 170, height = 85, units = "mm", dpi = 300)


# Connectivity and runout map ##################################################

(load("C:\\sda\\Workspace\\sedconnect\\runoutSim_wConnectFeature.Rd"))
library(ggnewscale) # Allow multiple scales in a map
library(ggspatial) # for scale bar

conn <- connToRaster(multi_sim_paths, dem)
vel <- velocityToRaster(multi_sim_paths, dem)
rel_paths <- rasterCdf(walksToRaster(multi_sim_paths, method = "freq", dem))


# Plot conn

df_hill <- as.data.frame(hill, xy = TRUE, na.rm = TRUE)
df_conn <- as.data.frame(conn, xy = TRUE, na.rm = TRUE)
df_path <- as.data.frame(rel_paths, xy = TRUE, na.rm = TRUE)
df_srcprob <- as.data.frame(src_pred, xy = TRUE, na.rm = TRUE)
df_srcarea <- as.data.frame(src_area, xy = TRUE, na.rm = TRUE)

m.srcprob <- ggplot() +
  geom_tile(data=df_hill, aes(x=x, y=y, fill = hillshade),
            show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +
  new_scale("fill") +
  geom_tile(data = df_srcprob, aes(x = x, y = y, fill = lyr1)) +
  scale_fill_viridis(name = "Source\nProbability\n", 
                     option = "viridis", direction = -1,
                     alpha = 0.7, na.value = "transparent") +
  geom_sf(data = drainage_network, 
          fill = alpha("#99d2ff", 0.5), 
          color = alpha("#99d2ff", 0.5), size = 0.8) +
  #geom_sf(data = bnd_catchment, fill = NA, color = "#566573" , size = 2)+
  xlab("") +
  ylab("") +
  coord_sf() +
  theme_light() +
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.5, "cm"),
        text = element_text(size = 9), 
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 6), 
        axis.text.y = element_text(angle = 90)) +
  annotation_scale(aes(style = "ticks", location = "br"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) 

m.srcarea <- ggplot() +
  geom_tile(data=df_hill, aes(x=x, y=y, fill = hillshade),
            show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +
  new_scale("fill") +
  geom_tile(data = df_srcarea, aes(x = x, y = y, fill = ""), alpha = 0.7) +
  scale_fill_manual(name = "Classified\nSource Area", values = "#e74c3c") +
  geom_sf(data = drainage_network, 
          fill = alpha("#99d2ff", 0.5), 
          color = alpha("#99d2ff", 0.5), size = 0.8) +
  #geom_sf(data = bnd_catchment, fill = NA, color = "#566573" , size = 2)+
  xlab("") +
  ylab("") +
  coord_sf() +
  theme_light() +
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.5, "cm"),
        text = element_text(size = 9), 
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 6), 
        axis.text.y = element_text(angle = 90)) +
  annotation_scale(aes(style = "ticks", location = "br"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) 

m.conn <- ggplot() +
  geom_tile(data=df_hill, aes(x=x, y=y, fill = hillshade),
            show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +
  new_scale("fill") +
  geom_tile(data = df_conn, aes(x = x, y = y, fill = connectivity_prob)) +
  scale_fill_viridis(name = "Connectivity\nProbability\n", 
                     option = "viridis", direction = -1,
                     alpha = 0.7, na.value = "transparent") +
  geom_sf(data = drainage_network, 
          fill = alpha("#99d2ff", 0.5), 
          color = alpha("#99d2ff", 0.5), size = 0.8) +
  #geom_sf(data = bnd_catchment, fill = NA, color = "#566573" , size = 2)+
  xlab("") +
  ylab("") +
  coord_sf() +
  theme_light() +
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.5, "cm"),
        text = element_text(size = 9), 
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 6), 
        axis.text.y = element_text(angle = 90)) +
  annotation_scale(aes(style = "ticks", location = "br"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) 


m.paths <- ggplot() +
  geom_tile(data=df_hill, aes(x=x, y=y, fill = hillshade),
              show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +
  new_scale("fill") +
  geom_tile(data = df_path, aes(x = x, y = y, fill = freq)) +
  scale_fill_viridis(name = "Traverse\nFrequency\n(Quantiles)", alpha = 0.7, 
                     direction = -1, na.value = "transparent") +
  geom_sf(data = drainage_network, 
          fill = alpha("#99d2ff", 0.5), 
          color = alpha("#99d2ff", 0.5), size = 0.8) +
  #geom_sf(data = bnd_catchment, fill = NA, color = "#566573" , size = 2)+
  xlab("") +
  ylab("") +
  coord_sf() +
  theme_light() +
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.5, "cm"),
        text = element_text(size = 9), 
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 6), 
        axis.text.y = element_text(angle = 90)) +
  annotation_scale(aes(style = "ticks", location = "br"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) 


m.paths_conn <- m.paths + m.conn

ggsave("Case_Study/Figures/connectivity_maps.png", plot = m.paths_conn, width = 170, height = 120, units = "mm", dpi = 300)

m.src <- m.srcprob + m.srcarea
ggsave("Case_Study/Figures/src_maps.png", plot = m.src, width = 170, height = 120, units = "mm", dpi = 300)


m.reg_results <- (m.srcprob + m.srcarea) / (m.paths + m.conn)

ggsave("Case_Study/Figures/regional_maps.png", plot = m.reg_results, width = 170, height = 240, units = "mm", dpi = 300)

