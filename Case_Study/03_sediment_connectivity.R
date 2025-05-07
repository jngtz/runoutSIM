library(runoutSim)
library(terra)
library(sf)


# Load data ####################################################################

# Load digital elevation model (DEM)
dem <- rast("Dev/Data/elev_nosinks.tif") # use sink filled DEM to remove pits and flats 

river_channel <- st_read("Dev/Data/river_channel.shp")
stream_channels <- st_read("Dev/Data/river_rio_olivares.shp")

bnd_catchment <- st_read("Dev/Data/basin_rio_olivares.shp")
# buffer stream channels
buffer_stream <- st_buffer(stream_channels, dist = 30)

# combine to one feature
drainage_network <- st_sf(st_union(st_union(st_geometry(river_channel), st_geometry(buffer_stream), is_coverage = TRUE)))

leafmap(drainage_network, color = "lightblue")

# rasterize drainage_network to view cells that will be considered for connectivity
r_dn = rasterize(drainage_network, dem)

source_grid <- rast("Data/auto_classified_source_areas.tif")

# View data
leafmap(bnd_catchment, 
        color = '#f7f9f9', 
        fill_color = '#FF000000', 
        weight = 4) %>%
  
  leafmap(source_grid, palette = list(classes = 1, colors = "pink", labels = "Source cells")) %>%
  
  leafmap(r_dn, 
          label = "Drainage Network", 
          palette = list(classes = 1, colors = "#99d2ff", labels = ""), 
          opacity = 0.8)


# Data pre-processing ##########################################################

# We create a connectivity feature to allow for quicker processing 
# it is essentially an index of all the cells covered by the feature max in 
# the given DEM.

feature_mask <- makeConnFeature(drainage_network, dem)

# Define source cells from source grid #########################################

sum(values(source_grid), na.rm = TRUE)

# ! make below into a function
# Find cells where the value is 1
source_cells <- which(values(source_grid) == 1)

# Extract the coordinates of these cells
source_xy <- xyFromCell(source_grid, source_cells)
source_xy <- source_xy[1:10,]

# Create a list of matricies for input to pcmRW
source_l <- list()
for(i in 1:nrow(source_xy)){
  source_l[[i]] <- matrix(source_xy[i,], ncol=2)
}

rw_l <- list()
for(i in 1:length(source_l)){
  print(i)
  rw_l[[i]] <- runoutSim(dem = dem, xy = source_l[[i]], mu = 0.08, md = 140, 
            slp_thresh = 50, exp_div = 3.0, per_fct = 1.95, walks = 1000,
            source_connect = TRUE, connect_feature = feature_mask)
}


trav_freq <- walksToRaster(rw_l, dem)
trav_prob <- rasterCdf(trav_freq)
