# LIBRARIES ####################################################################
setwd("")
setwd("C:\\Users\\ku52jek\\Documents\\GitHub\\sedconnect\\R")
source("run_connect.R")
source("~/R/sedconnect/R/run_connect.R")


library(rgdal)
library(raster)
library(terra)
library(doParallel)
library(foreach)
library(rgeos)
library(runoptGPP)

# TASKS
#  x add crop to runConnect to reduce size of DEM: control with extent 
#  o Make blockElvCropExtent faster using foreach
#  x add path_list[[1]] = cntr_cell to rw_pcm.R and random_walk.R

# NOTES
# - terra() functions are faster than R, but cannot be serialized (shared across nodes)
#   one way around this is create the terra objects withn the foreach loop.

# Applications
# Determine avalanche source areas that could impact forest, ski-resorts etc...


# FUNCTIONS ####################################################################



blockCropExtent <- function(block_xy, dem, extend_crop = 6000){
  # Determine Extent for fast processing
  # Requires UTM/local projection to perform
  # extend_crop can be set to max. expected runout
  
  sp_xys <- sapply(block_xy, FUN = SpatialPoints, proj4string = crs(dem))
  
  extents <- block_xy
  
  for(i in 1:length(block_xy)){
    
    for(k in 1:length(sp_xys[[i]])){
      extents[[i]] <- list()
    }
  }
  
  # Calc crop extent of each point using max elev
  for(i in 1:length(block_xy)){
    
    for(k in 1:length(sp_xys[[i]])){
      extents[[i]][[k]] <- extent(sp_xys[[i]][k]) + extend_crop # calc extent of each point
    }
  }
  
  return(extents)
  
}

blockElvCropExtent <- function(block_xy, block_cell, dem, extend_crop = 6000){
  # Determine Extent for fast processing
  # Requires UTM/local projection to perform
  # extend_crop can be set to max. expected runout
  # o make this parallel to run faster
  
  max_elvs <- sapply(block_cell, FUN = extract, x = dem )
  sp_xys <- sapply(block_xy, FUN = SpatialPoints, proj4string = crs(dem))
  
  extents <- block_xy
  
  for(i in 1:length(block_xy)){
    
    for(k in 1:length(sp_xys[[i]])){
      extents[[i]] <- list()
    }
  }
  
  # Calc crop extent of each point using max elev
  for(i in 1:length(block_xy)){
    
    for(k in 1:length(sp_xys[[i]])){
      #extents[[i]][[k]] <- extent(sp_xys[[i]][k]) # calc extent of each point
      
      #dem_crop <- crop(dem, extent(sp_xys[[i]][k]) + extend_crop)
      #dem_crop[dem_crop > max_elvs[[i]][k] + 10] <- NA # + 10 m in elevation
      #dem_trim <- trim(dem_crop)
      #dem_extent <- extent(dem_trim)
      
      dem_crop <- terra::crop(test, extent(sp_xys[[i]][k]) + extend_crop)
      dem_crop[dem_crop > max_elvs[[i]][k] + 10] <- NA
      dem_crop <- terra::trim(dem_crop)
      
      
      #dem_extent <- extent(raster(dem_crop))
      
      
      dem_extent <- ext(dem_crop)
      
      #test_crop <- crop(dem_crop, dem_extent)
      #plot(test_crop)
      #plot(sp_xys[[i]][k], add = TRUE)
      
      extents[[i]][[k]] <- extent(raster(dem_extent))
      
    }
  }
  
  return(extents)
  
}

runConnectCrop <- function(crop_extent,  dem, object, xy, mu = 0.1, md = 40, slp_thresh = 30, exp_div = 3, per_fct = 2, reps = 100){
  # For faster processing on large DEMs.
  
  #st <- Sys.time()
  dem_crop <- terra::crop(rast(dem), crop_extent)
  object_crop <- terra::crop(rast(object), crop_extent) 
  
  #dem_crop <- raster(dem_crop)
  #object_crop <- raster(object_crop)
  
  if(is.na(minValue(object_crop))){
    # assign probability value of 0 if no object in extent
    prob <- 0
    
  } else {
    
    prob <- runConnect(dem_crop, object_crop, xy, mu , md , slp_thresh , exp_div, per_fct , reps)
    
  }
  #Sys.time() - st
  
  return(prob)
  
}


blockTargetXYs <- function(source_grid, target_cells){
  sapply(target_cells, FUN = xyFromCell, object = source_grid)
}

blockTargetCells <- function(source_grid, n_cores = 10){
  # Create blocks of cells for parallel processing
  
  target_cells = which(values(source_grid)==1)
  n_target_cells = length(target_cells)
  
  cells_split <- round(n_target_cells/n_cores)
  
  group_cells <- list()
  group_cells[[1]] <- c(1, cells_split)
  
  for(i in 1:(n_cores)){
    start_cell <- (i-1)*cells_split+1
    end_cell <- (i)*cells_split
    group_cells[[i]] <-c(start_cell, end_cell)
  }
  
  group_cells[[n_cores]] <- c(group_cells[[n_cores]][1], n_target_cells)
  
  cell_blocks <- list()
  
  for(i in 1:(n_cores)){
    start_cell <- group_cells[[i]][1]
    end_cell <- group_cells[[i]][2]
    cell_blocks[[i]] <- target_cells[start_cell:end_cell]
  }
  
  return(cell_blocks)
}


blockExtentRunConnect <- function(block, block_xy, block_extent,
                                  dem, object, mu = 0.1, md = 40, slp_thresh = 30, exp_div = 3, per_fct = 2, reps = 100){
  
  block_result <- rep(NA, times = nrow(block_xy[[block]]))
  
  for(i in 1:length(block_result)){
    
    xy <- block_xy[[block]][i,]
    block_result[i] <- runConnectCrop(block_extent[[block]][[i]], dem, object, xy, mu, md, slp_thresh, exp_div, per_fct, reps)
    
  }
  
  return(block_result)
  
} #block function


blockRunConnect <- function(block, block_xy,
                            dem, object, mu = 0.1, md = 40, slp_thresh = 30, exp_div = 3, per_fct = 2, reps = 100){
  
  block_result <- rep(NA, times = nrow(block_xy[[block]]))
  
  for(i in 1:length(block_result)){
    
    xy <- block_xy[[block]][i,]
    block_result[i] <- runConnect(dem, object, xy, mu, md, slp_thresh, exp_div, per_fct, reps)
    
  }
  
  return(block_result)
  
} #block function


# PARALLEL PROCESSING PARAMETERS ################################################

# Number of cores available ('logical processors')
n_cores <- 8

# LOAD RASTER DATA #############################################################

#setwd("~/Data/m_buchhart")
setwd("~/Desktop/sda/sedconnect/data")
# Read raster using raster package
dem <- raster("elev.tif")

# Read Polygons and Source Points using RGDAL package
runout_polygons <- readOGR(".", "debris_flow_runout_polygons")
runout_polygons$objectid <- 1:length(runout_polygons)

source_points <- readOGR(".", "debris_flow_source_points")

# Read river polyline
river_channel <- readOGR(".", "river_channel")
river <- readOGR(".", "river_rio_olivares")
river_buff <- buffer(river, 30)
object = rasterize(river_buff, dem)


# THRESHOLD SOURCE PREDICTION ##################################################

source_grid<- raster("auto_classified_source_areas.tif") # from "02_source_area_prediction.R"

# # For quick testing #
# ext_test <- extent(source_grid)
# ext_test@ymax = 6327000
# ext_test@ymin = 6320000
# ext_test@xmin = 387000
# ext_test@xmax = 396000
# 
# dem_crop <- crop(dem, ext_test)
# source_grid_crop <- crop(source_grid, ext_test)
# object_crop <- crop(object, ext_test)
# 
# plot(dem_crop)
# plot(source_grid_crop)
# plot(object_crop)
# 
# dem <- dem_crop
# source_grid <- source_grid_crop
# object <- object_crop

# RUN GPP PER CELL FOR SUBCATCHMENT ############################################

# Determine target cell locations
block_cell <- blockTargetCells(source_grid, n_cores)
block_xy <- blockTargetXYs(source_grid, block_cell)
#block_extent <- blockCropExtent(block_xy, dem, extend_crop = 6000)

# Visualize blocks
plot(dem)
for(i in 1:n_cores)plot(SpatialPoints(block_xy[[i]]), add = T, col = i)

# Run GPP for each target cell in parallel using FOREACH
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)


# Turn below itno
connect_probs <-  foreach(BLOCK=1:n_cores, .combine=c,
                          .packages=c('rgdal','raster', 'terra', 'rgeos')) %dopar% {
                            
                            blockRunConnect(block = BLOCK, block_xy, 
                                            #block_extent,
                                            dem, object, 
                                            mu = 0.08, 
                                            md = 140, 
                                            slp_thresh = 35, 
                                            exp_div = 3, 
                                            per_fct = 1.95, 
                                            reps = 1000)
                          }


parallel::stopCluster(cl)


save(connect_probs, file = "rio_olivares_srcconn_probs_vector.Rd")
prob_grid <- source_grid
values(prob_grid) <- NA
target_cells <- unlist(block_cell)

prob_grid[target_cells] = rnorm(length(target_cells)) # test with random grid
prob_grid[target_cells] = connect_probs

writeRaster(prob_grid, filename = "rio_olivares_srcconn_probs.tif", format = "GTiff", overwrite = TRUE )


# Make figures ###############################################################

library(ggplot2) # Plotting maps
library(patchwork) # Tile multiple maps together
library(ggnewscale) #
library(ggspatial)
library(sf)

# Functions for pretty mapping

rasterToDF <- function(x){
  # creates a data frame from a raster object
  # helps for mapping with ggplot2
  df <- as.data.frame(x, xy = TRUE)
  df <- df[!is.na(df[,3]),]
  names(df) <- c("x", "y", "z")
  df
}


ggMultiHillshade <- function(slope, aspect, angles = c(70, 60, 55), directions = c(350, 15, 270),
                             alphas = c(0.7, 0.5, 0.65)){
  # Creates multiple hillshade object for ggplot
  
  hs_top <-  hillShade(slope, aspect, angle=angles[1], direction=directions[1])
  hs_middle <-  hillShade(slope, aspect, angle=angles[2], direction=directions[2])
  hs_bottom <-  hillShade(slope, aspect, angle=angles[3], direction=directions[3])
  
  
  hs_top_df <- rasterToDF(hs_top)
  hs_middle_df <- rasterToDF(hs_middle)
  hs_bottom_df <- rasterToDF(hs_bottom)
  
  gg_hillshades <-  ggplot() +
    geom_sf() +
    geom_raster(data=hs_bottom_df, aes(x=x, y=y, fill = z), show.legend = FALSE, alpha = alphas[1]) +
    scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +
    
    ggnewscale::new_scale("fill") +
    geom_raster(data=hs_middle_df, aes(x=x, y=y, fill = z), show.legend = FALSE, alpha = angles[2]) +
    scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +
    
    ggnewscale::new_scale("fill") +
    geom_raster(data=hs_top_df, aes(x=x, y=y, fill = z), show.legend = FALSE, alpha = angles[3]) +
    scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF")
  
  return(gg_hillshades)
  
}

geom_hillshades <- ggMultiHillshade(slope, aspect, angles = c(70, 60, 55), directions = c(350, 15, 270),
                                    alphas = c(0.7, 0.5, 0.65))



# Read river polyline
river_channel <- readOGR(".", "river_channel")
river <- readOGR(".", "river_rio_olivares")
river_buff <- buffer(river, 30)
object = rasterize(river_buff, dem)

rivers <- st_read("river_rio_olivares.shp")
runout <- st_read("debris_flow_runout_polygons.shp")
bnd <- st_read("basin_rio_olivares.shp")
#Clip rivers
rivers <- st_intersection(rivers, bnd )


parea <- raster("gpp_process_area.tif") # from "03b_pcm_runout_optimization.R"
parea_cdf <- rasterCdf(parea)

conn <- raster("src_connect_max.tif") #change to "rio_olivares_srcconn_probs.tif" "max" just for illustration


# For visualization make a hillshade model from the DEM
slope <- terrain(dem, opt='slope')
aspect <- terrain(dem, opt='aspect')
hillshade <- hillShade(slope, aspect, angle=40, direction=270)

# Transform raster data into a data frame for use with ggplot
hillshade_df <- as.data.frame(hillshade, xy = TRUE)
hillshade_df <- hillshade_df[!is.na(hillshade_df[,3]),]


gppGGplot <- function(x){
  #parea_cdf <- rasterCdf(x)
  x <- as.data.frame(x, xy = TRUE)
  x <- x[!is.na(x[,3]),]
  return(x)
}

parea_df<- gppGGplot(parea_cdf)
conn_df<- gppGGplot(conn)


map.parea <- geom_hillshades +
  
  geom_sf(data = bnd, alpha = 0.4, fill = "white") +
  
  new_scale("fill") +
  geom_tile(data=parea_df, aes(x=x, y=y, fill = parea_df[,3]) ) +
  scale_fill_viridis_c(name = "Runout\nfreq.\n(quantiles)",  alpha = 0.6, direction = -1) +
  
  geom_sf(data = rivers, colour = "#85C1E9", fill = alpha("#85C1E9", 0.5), lwd = 0.3) +
  #geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  geom_sf(data = bnd, colour = "#7f8c8d", fill = NA) +
  
  #geom_sf(data = runout, color=alpha("black", 0.3), fill = alpha("white",0.2) , size = 0.3) +
  
  xlab("") +
  ylab("") +
  theme_bw() +
  annotation_scale(aes(style = "ticks", location = "br"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) +
  theme(text = element_text(size = 7), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6), axis.text.y = element_text(angle = 90),
        legend.key.width=unit(0.5,"cm"), legend.key.height = unit(0.5, "cm"))



map.conn <- geom_hillshades +
  
  geom_sf(data = bnd, alpha = 0.4, fill = "white") +
  
  new_scale("fill") +
  geom_tile(data=conn_df, aes(x=x, y=y, fill = conn_df[,3]) ) +
  scale_fill_viridis_c(name = "Source\nconnectivity\n(prob.)",  alpha = 0.6, direction = -1) +
  
  geom_sf(data = rivers, colour = "#85C1E9", fill = alpha("#85C1E9", 0.1), lwd = 0.3) +
  #geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  geom_sf(data = bnd, colour = "#7f8c8d", fill = NA) +
  
  #geom_sf(data = runout, color=alpha("black", 0.3), fill = alpha("white",0.2) , size = 0.3) +
  
  xlab("") +
  ylab("") +
  theme_bw() +
  annotation_scale(aes(style = "ticks", location = "br"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) +
  theme(text = element_text(size = 7), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6), axis.text.y = element_text(angle = 90),
        legend.key.width=unit(0.5,"cm"), legend.key.height = unit(0.5, "cm"))

#setwd("C:/Users/ku52jek/Nextcloud/project/results/figures/report")

map.parea + ggtitle("Simulated runout") +
  map.conn + ggtitle("Source connectivity") + 
  plot_annotation(tag_levels = "A")
ggsave("map__olivares_runfreq_src_conn.png", dpi = 300, width = 7.4, height = 5, units = "in")


map.conn
ggsave("map__olivares_srcconn_index.png", dpi = 300, width = 5, height = 6, units = "in")


# FOR EGU #######################################################################

sub_catchments <- readOGR("sub_catchments.shp")
sub_catchment <- sub_catchments[sub_catchments$ID == 56,]

bnd_sub <- st_as_sf(sub_catchment)
river_sub <- st_intersection(bnd_sub, rivers)

map.conn_egu <- ggplot() +
  geom_sf() +
  geom_raster(data=hillshade_df, aes(x=x, y=y, fill = layer),
              show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +
  
  geom_sf(data = bnd_sub, alpha = 0.4, fill = "white") +
  
  new_scale("fill") +
  geom_tile(data=conn_df, aes(x=x, y=y, fill = conn_df[,3]),
              show.legend = TRUE) +
  scale_fill_viridis_c(name = "Source\nconnectivity\n(prob.)",  alpha = 0.6, direction = -1) +
  
  geom_sf(data = river_sub, colour = "#85C1E9") +
  #geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  geom_sf(data = bnd_sub, colour = "#7f8c8d", fill = NA) +
  
  annotation_scale(aes(style = "ticks", location = "br"), text_family = "Arial", text_cex = 0.6,
                   bar_cols = "black", line_width = 0.7) +
  
  xlab("") +
  ylab("") +
  theme_void() +
  theme(text = element_text(family = "Arial", size = 7),
        axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "right")

map.conn_egu
ggsave("map_egu_olivares_srcconn_index.png", dpi = 300, width = 5, height = 5.5, units = "in")

# PRETTY FOR EGU #########################




library(metR)

map.conn.egu <- geom_hillshades +
  
  geom_sf(data = bnd, alpha = 0.4, fill = "white") +
  
  new_scale("fill") +
  geom_tile(data=conn_df, aes(x=x, y=y, fill = conn_df[,3]) ) +
  scale_fill_viridis_c(name = "Source\nconnectivity\n(prob.)",  alpha = 0.6, direction = -1) +
  
  geom_sf(data = rivers, colour = "#85C1E9", fill = alpha("#85C1E9", 0.1), lwd = 0.3) +
  #geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  geom_sf(data = bnd, colour = "#7f8c8d", fill = NA) +
  
  #geom_sf(data = runout, color=alpha("black", 0.3), fill = alpha("white",0.2) , size = 0.3) +
  
  xlab("") +
  ylab("") +
  theme_bw() +
  annotation_scale(aes(style = "ticks", location = "br"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) +
  theme(text = element_text(size = 7), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6), axis.text.y = element_text(angle = 90),
        legend.key.width=unit(0.5,"cm"), legend.key.height = unit(0.5, "cm"))


map.conn.egu
ggsave("map_egu_olivares_srcconn_index.png", dpi = 300, width = 5, height = 5.5, units = "in")

# Export runout and connectivity as KML for Google Earth ################################

library(viridisLite)
p_conn <- projectRaster(conn, crs="+proj=longlat +datum=WGS84", method='ngb')
KML(p_conn, filename = "srconn_.kml", col = viridis(256, direction = -1), maxpixels = ncell(conn), overwrite = TRUE)

p_parea <- projectRaster(parea_cdf, crs="+proj=longlat +datum=WGS84", method='ngb')
KML(p_parea, filename = "parea_cdf_.kml", col = viridis(256, direction = -1), maxpixels = ncell(conn), overwrite = TRUE)



