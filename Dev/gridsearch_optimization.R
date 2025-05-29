# Development (Dev) environment for optimizing random walk runout simulations
library(runoutSim)
library(runoptGPP)

# New strategy
# - 1st estimate runout distance PCM (quicker)
# - 2nd estimate runout path (shorter runout will help reduce processing times.)


# If having bounding box length issue ... adjust pcmPerformance to also calculate
# length as the distance from the highest elevation point (source point)
# to the distance of the lowest elevation points of the runout.

# OR calculate D8 distance, the distance along the steepest path


# Optimizing an individual runout path simulation ##############################

library(raster)
library(terra)
library(sf)

# Load digital elevation model (DEM)
dem <- raster("Dev/Data/elev_fillsinks_WangLiu.tif")
#dem <- raster("Dev/Data/elev_nosinks.tif") # use sink filled DEM to remove pits and flats 

# Load runout source points and polygons
source_points <- st_read("Dev/Data/debris_flow_source_points.shp")
source_points$run_id <- 1:nrow(source_points)
runout_polygons <- st_make_valid(st_read("Dev/Data/debris_flow_runout_polygons.shp"))
# ^ Need to clean this up so make valid not needed

# Select a single debris flow and source point for the example
#runout_polygons <- runout_polygons[10:20,]
runout_polygons$run_id <- 1:nrow(runout_polygons)

# Ad spatial join to associate runout with source with points

rungeom <- runoptGPP::runoutGeom(as_Spatial(runout_polygons), dem)
rungeom$fid <- rungeom$id <- NULL

runout_polygons <- cbind(runout_polygons, rungeom)

# Get corresponding source point
source_points  <- st_filter(st_as_sf(source_points), st_as_sf(runout_polygons))

leafmap(runout_polygons) %>% leafmap(source_points, col = "red")


## Set-up working environment to allow for parallel processing #################
library(doParallel)
library(foreach)

# DEFINE RW AND PCM GRID SEARCH SPACE ##########################################


pcmmu_vec <- seq(0.04, 0.6, by=0.01)
polyid_vec <- 1:nrow(source_points)


# PCM GRIDSEARCH OPTIMIZATION W PARALLELIZATION  ################################
setwd("C:\\sda\\Workspace\\sedconnect\\backup")

n_cores <- detectCores() -2
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

pcm_gridsearch_multi <-
  foreach(poly_id=polyid_vec, .packages=c('terra','raster', 'ROCR', 'sf', 'runoptGPP', 'runoutSim')) %dopar% {
    
    pcmGridsearch(dem,
                  slide_plys = runout_polygons, slide_src = source_points, slide_id = poly_id,
                  rw_slp = 40, rw_ex = 3, rw_per = 1.9, # from Goetz et al 2021
                  pcm_mu_v = pcmmu_vec, pcm_md_v = 40,
                  gpp_iter = 1000,
                  buffer_ext = NULL, buffer_source = 20,
                  predict_threshold = 0.5, save_res = TRUE,
                  plot_eval = FALSE,
                  saga_lib = NULL)
    
  }

parallel::stopCluster(cl)

# Get PCM optimal parameter set

pcm_gridsearch_multi <- list()
files <- list.files(pattern = "result_pcm_gridsearch_")
for (i in 1:length(files)) {
  res_nm <- paste("result_pcm_gridsearch_", i, ".Rd", 
                  sep = "")
  res_obj_nm <- load(res_nm)
  result_pcm <- get(res_obj_nm)
  pcm_gridsearch_multi[[i]] <- result_pcm
}


pcm_opt <- pcmGetOpt(pcm_gridsearch_multi, performance = "relerr", measure = "median", plot_opt = TRUE, from_save = TRUE)

# Create a plot of median relerr for mu
pcm_opt

grid_res <- data.frame(
  mu = as.numeric(rownames(pcmGetGrid(performance = "relerr", measure = "median", from_save = TRUE))),
  relerr_median = pcmGetGrid(performance = "relerr", measure = "median", from_save = TRUE)[,1],
  error_median = pcmGetGrid(performance = "error", measure = "median", from_save = TRUE)[,1],
  abs_error_median = abs(pcmGetGrid(performance = "error", measure = "median", from_save = TRUE)[,1]),
  iqr = pcmGetGrid(performance = "relerr", measure = "IQR", from_save = TRUE)[,1]
)

library(ggplot2)
library(dplyr)


grid_res <- grid_res %>%
  mutate(
    ymin = relerr_median - iqr / 2,
    ymax = relerr_median + iqr / 2
  )

# Create the plot with IQR as ribbon
ggplot(grid_res, aes(x = mu, y = relerr_median)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "gray80", alpha = 0.5) +
  geom_line(color = "#19191a", linewidth = 0.9) +
  geom_point(color = "#19191a", size = 1.5) +
  scale_x_continuous(breaks = seq(0, 2, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(
    x = expression("Sliding friction coefficient (" * mu * ")"),
    y = "Median relative runout distance error"
  ) +
  theme_classic(base_size = 10)



save(pcm_opt, file = "pcm_opt_params.Rd")
save(pcm_gridsearch_multi, file = "pcm_gridsearch_multi.Rd")

# PCM PARAM VALIDATION W SPATIAL CV #############################################

pcm_spcv <- pcmSPCV(pcm_gridsearch_multi, slide_plys = runout_polygons,
                    n_folds = 5, repetitions = 10, from_save = FALSE)

freq_pcm <- pcmPoolSPCV(pcm_spcv, plot_freq = TRUE)
freq_pcm

# Save results
setwd("C:\\sda\\Workspace\\sedconnect")
save(global_run, file = "gridsearch_optimization_RunoutSim.Rd")
(load("gridsearch_optimization_RunoutSim.Rd"))



# Retrieve overall performance

results_df <- foreach(
  i = seq_len(nrow(runout_polygons)),
  .packages=c('terra','raster', 'ROCR', 'sf', 'runoptGPP', 'runoutSim')) %dopar% {
    out <- 
      pcmPerformance(
        dem            = dem,
        slide_plys     = runout_polygons,
        slide_src      = source_points,
        slide_id       = i,
        rw_slp         = 40,
        rw_ex          = 3,
        rw_per         = 1.9,
        pcm_mu         = 0.08,
        pcm_md         = 40,
        gpp_iter       = 1000,
        buffer_ext     = NULL,
        buffer_source  = 20,
        plot_eval      = FALSE,
        return_features= TRUE
      )
  }

results_df

stopCluster(cl)


# Extract error
runout_polygons$sim_roc <- sapply(results_df, function(x) x$roc)
runout_polygons$sim_length_relerror <- sapply(results_df, function(x) x$length.relerr)
runout_polygons$sim_length_error <- sapply(results_df, function(x) x$length.error)
runout_polygons$sim_length <- runout_polygons$length + runout_polygons$sim_length_error 

median(runout_polygons$sim_length_relerror)
IQR(runout_polygons$sim_length_relerror)


median(runout_polygons$sim_length_error)


png(filename="hist_error_plots.png", res = 300, width = 7.5, height = 6,
    units = "in", pointsize = 11)

par(family = "Arial", mfrow = c(2,2), mar = c(4, 3, 0.5, 0.5),
    mgp = c(2, 0.75, 0))

hist(runout_polygons$sim_roc, breaks = 20,
     main = "",
     xlab = "Runout path AUROC")

hist(runout_polygons$sim_length_relerror, breaks = 40,,
     main = "",
     xlab = "Runout length relative error")

hist(runout_polygons$sim_length_error, breaks = 40,,
     main = "",
     xlab = "Runout length error (m)")

plot(runout_polygons$sim_length, runout_polygons$length,
     xlab = "Simulated runout length (m)",
     ylab = "Observed runout length (m)", pch = 20)

dev.off()



# > change below... 

# Map results ##################################################################












# Get buffered grid cells used in model training

sp_runout_polygons <- as_Spatial(runout_polygons)
sp_runout_polygons$rowid <- 1:nrow(sp_runout_polygons)

sp_source_points <- as_Spatial(source_points)

mbo_source_xy <- list()
mbo_source_grid <- list()
for(i in 1:nrow(runout_polygons)){
  
  runout_single <- sp_runout_polygons[i,]
  
  sel_over_start_point  <- sp::over(sp_source_points, runout_single)
  sel_start_point <- sp_source_points[!is.na(sel_over_start_point$rowid),]
  
  source_buffer <- sf::st_buffer(st_as_sf(sel_start_point), dist = 20)
  source_grid <- raster::rasterize(source_buffer, dem, field=1 )
  source_grid <- raster::mask(source_grid, runout_single )
  
  source_cells <- which(values(source_grid) == 1)
  
  # Extract the coordinates of these cells
  mbo_source_xy[[i]] <- xyFromCell(source_grid, source_cells)
  mbo_source_grid[[i]] <- rast(source_grid)
}

mbo_source_grid <- rast(mbo_source_grid)
mbo_source_grid <- app(mbo_source_grid, fun = sum, na.rm = TRUE)
# Find cells where the value is 1
source_cells <- which(values(mbo_source_grid) == 1)

# Extract the coordinates of these cells
source_xy <- xyFromCell(dem, source_cells)

source_l <- list()
for(i in 1:nrow(source_xy)){
  source_l[[i]] <- matrix(source_xy[i,], ncol=2)
}


library(parallel)
# Define number of cores to use
n_cores <- detectCores() -2

dem_terra <- terra::rast('C:\\sda\\GitProjects\\runoutSim\\Dev/Data/elev_fillsinks_WangLiu.tif')
#dem_terra <- terra::rast('C:\\sda\\GitProjects\\runoutSim\\Dev\\Data\\elev.tif')


packed_dem <- wrap(dem_terra)

# Create parallel loop
cl <- makeCluster(n_cores, type = "PSOCK") # Open clusters

# Load objects and 'custom' functions to each cluster
# clusterExport(cl, varlist = c("runoutSim", "euclideanDistance", "adjCells",
#                              "adjRowCol", "pcm", "packed_dem","sourceConnect",
#                              "feature_mask"))

clusterExport(cl, varlist = c("packed_dem"))

# Load required packages to each cluster
clusterEvalQ(cl, {
  library(terra)
  library(runoutSim)
})


# Load objects and 'custom' functions to each cluster
# clusterExport(cl, varlist = c("runoutSim", "euclideanDistance", "adjCells",
#                              "adjRowCol", "pcm", "packed_dem","sourceConnect",
#                              "feature_mask"))





multi_sim_paths <- parLapply(cl, source_l, function(x) {
  
  runoutSim(dem = unwrap(packed_dem), xy = x, 
            mu = 0.08, 
            md = 40, 
            slp_thresh = 40, 
            exp_div = 3, 
            per_fct = 1.9, 
            walks = 1000)
})



stopCluster(cl) 

trav_freq <- walksToRaster(multi_sim_paths, dem_terra)
trav_prob <- runoutSim::rasterCdf(trav_freq)

trav_vel <- velocityToRaster(multi_sim_paths, dem_terra)

leafmap(runout_polygons) %>%
  leafmap(trav_freq) %>%
  leafmap(trav_prob) %>%
  leafmap(trav_vel, palette = 'plasma') %>%
  leafmap(source_points, color = "red") 

########################################################



# Extract all gpp.parea rasters
raster_list <- lapply(results_df, function(x) rast(x$gpp.parea))

# Step 1: Merge to find the full union extent
target_ext    <- ext(dem)

# Step 2: Extend each raster to the same extent (fills gaps with NA)
aligned_list <- lapply(raster_list, function(r) {
  extend(r, target_ext)
})

# Step 3: Stack and sum
aligned_stack <- rast(aligned_list)

# Step 4: Add pixel values, ignoring NAs
trav_freq <- app(aligned_stack, fun = sum, na.rm = TRUE)

trav_prob <- runoutSim::rasterCdf(trav_freq)


leafmap(trav_prob) %>% leafmap(runout_polygons) %>% leafmap(source_points, color = "red")


plot(global_run)