# runoutSIM case study: Runout optimization ####################################

# This script is 1 of 3 scripts presenting a case study (demo) for the runoutSIM
# package for regional mass movement runout simulation.

# It provides an example of using runoutSIM for regional debris flow runout
# model optimization using data for the Río Olivares basin of the semi-arid 
# central Chilean Andes. 

# By Jason Goetz, PhD (jgoetz@wlu.ca)
# Department of Geography and Environmental Studies,
# Wilfrid Laurier University, Canada

# August 7, 2025

## Steps: ####
# 1. Load required packages
# 2. Load spatial data
# 3. Perform runout model optimization
# 4. Get optimal runout parameters
# 5. Validate model performance
# 6. Run model from training source
# 7. Coerce modelling results to rasters (export)

# Load required packages #######################################################

# For runout path simulation and optimization
library(runoutSim)
#remotes::install_github("jngtz/runoptGPP")
library(runoptGPP)

# For spatial data handling
library(raster)
library(terra)
library(sf)

# For parallel processing
library(doParallel)
library(foreach)
library(parallel)
library(future.apply)


# For data handling and visualization
library(ggplot2)
library(dplyr)

# For random grid search optimization
library(lhs)


# Load spatial data  ###########################################################

# Load digital elevation model (DEM)
dem <- raster("Data/elev_fillsinks_WangLiu.tif")

# Hillshade for visualization 
slope <- terrain(rast(dem), "slope", unit="radians")
aspect <- terrain(rast(dem), "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270)


# Load runout source points and polygons
source_points <- st_read("Data/debris_flow_source_points.shp")
source_points$run_id <- 1:nrow(source_points)
runout_polygons <- st_read("Dev/Data/debris_flow_runout_polygons.shp")
runout_polygons$run_id <- 1:nrow(runout_polygons)

# Add spatial join to associate run out with source with points

# Calculate runout polygon geometry
rungeom <- round(runoptGPP::runoutGeom(as_Spatial(runout_polygons), dem),2)
rungeom$fid <- rungeom$id <- NULL

runout_polygons <- cbind(runout_polygons, rungeom)

# Get corresponding source point
source_points  <- st_filter(st_as_sf(source_points), st_as_sf(runout_polygons))

## Visualize input data ####

leafmap(runout_polygons, opacity = 0.6) %>%
  leafmap(source_points, color = '#e74c3c') %>%
  leafmap(hill, palette = grey(0:100/100), opacity = 1, add_legend = FALSE) %>%
  leafmap(dem, palette = viridis::mako(10), opacity = 0.6)

# Perform runout model optimization ############################################

## Parallel setup ####

plan(multisession, workers = parallel::detectCores() - 2)

## Define grid search space ####
n_samples <- 50

# Generate Latin Hypercube Samples for 5 parameters: 3 RW + 2 PCM
lhs_combined <- randomLHS(n_samples, 5)

combined_samples <- data.frame(
  # RW parameters
  rw_slp = qunif(lhs_combined[,1], 20, 40),
  rw_ex  = qunif(lhs_combined[,2], 1.3, 3),
  rw_per = qunif(lhs_combined[,3], 1.5, 2),
  
  # PCM parameters
  pcm_mu = qunif(lhs_combined[,4], 0.05, 0.4),
  pcm_md = qunif(lhs_combined[,5], 20, 120)
)

save(combined_samples, file = "Dev/random_search_latin_hypercube.Rd")
## Build an evaluation function ####

evaluate_combined <- function(rw, pcm) {
  n_slides <- nrow(runout_polygons)
  
  # Pre-allocate results matrix
  results <- matrix(NA, nrow = 3, ncol = n_slides,
                    dimnames = list(c("roc", "length_error", "length_relerr"), NULL))
  
  for (i in 1:n_slides) {
    out <- tryCatch({
      pcmPerformance(
        dem            = dem,
        slide_plys     = runout_polygons,
        slide_src      = source_points,
        slide_id       = i,
        rw_slp         = rw$rw_slp,
        rw_ex          = rw$rw_ex,
        rw_per         = rw$rw_per,
        pcm_mu         = pcm$pcm_mu,
        pcm_md         = pcm$pcm_md,
        gpp_iter       = 1000,
        buffer_source  = NULL,
        plot_eval      = FALSE,
        return_features= FALSE
      )
    }, error = function(e) NULL)
    
    if (is.null(out)) {
      results[, i] <- c(roc = 0, length_error = Inf, length_relerr = Inf)
    } else {
      results[, i] <- c(roc = out$roc, 
                        length_error = out$length.error, 
                        length_relerr = out$length.relerr)
    }
  }
  
  # Compute mean and median over all slides
  mean_roc <- mean(results["roc", ], na.rm = TRUE)
  mean_len <- mean(results["length_error", ], na.rm = TRUE)
  mean_relerr <- mean(results["length_relerr", ], na.rm = TRUE)
  
  median_roc <- median(results["roc", ], na.rm = TRUE)
  median_len <- median(results["length_error", ], na.rm = TRUE)
  median_relerr <- median(results["length_relerr", ], na.rm = TRUE)
  
  return(list(
    mean_roc = mean_roc,
    mean_length_error = mean_len,
    mean_relerr = mean_relerr,
    median_roc = median_roc,
    median_length_error = median_len,
    median_relerr = median_relerr
  ))
}

## Perform / run random search ####

combined_results <- future_lapply(1:nrow(combined_samples), function(k) {
  row <- as.list(combined_samples[k, ])
  
  res <- evaluate_combined(row, row)  # rw and pcm are from the same row
  
  c(sample_idx = k,
    rw_slp = row$rw_slp,
    rw_ex  = row$rw_ex,
    rw_per = row$rw_per,
    pcm_mu = row$pcm_mu,
    pcm_md = row$pcm_md,
    mean_roc = res$mean_roc,
    mean_length_error = res$mean_length_error,
    mean_relerr = res$mean_relerr,
    median_roc = res$median_roc,
    median_length_error = res$median_length_error,
    median_relerr = res$median_relerr)
})

# Close clusters
future:::ClusterRegistry("stop")

# Convert results to a data frame
combined_df <- do.call(rbind, lapply(combined_results, function(x) as.data.frame(as.list(x))))


# --- Optional: compute a weighted score ---
weight_len <- 2
combined_df$score <- combined_df$mean_roc - weight_len * combined_df$mean_length_error / max(combined_df$mean_length_error, 1e-6)

# --- Find best combination ---
best_row <- combined_df[which.max(combined_df$score), ]
best_params <- combined_df[which.max(combined_df$score), c("rw_slp","rw_ex","rw_per","pcm_mu","pcm_md")]

save(combined_df, file = "Dev/random_grid_search.Rd")

(load("Dev/random_grid_search.Rd"))

# Assess overall performance ###################################################

# Retrieve overall performance

n_cores <- detectCores() -2
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

results_df <- foreach(
  i = seq_len(nrow(runout_polygons)),
  .packages=c('terra','raster', 'ROCR', 'sf', 'runoptGPP', 'runoutSim'),
  .export = c("dem", "runout_polygons", "source_points", "best_params")) %dopar% {
    out <- 
      pcmPerformance(
        dem            = dem,
        slide_plys     = runout_polygons,
        slide_src      = source_points,
        slide_id       = i,
        rw_slp         = best_params$rw_slp,
        rw_ex          = best_params$rw_ex,
        rw_per         = best_params$rw_per,
        pcm_mu         = best_params$pcm_mu,
        pcm_md         = best_params$pcm_md,
        gpp_iter       = 1000,
        buffer_ext     = NULL,
        buffer_source  = 20,
        plot_eval      = FALSE,
        return_features= TRUE
      )
  }

results_df

parallel::stopCluster(cl)

save(results_df, file = "Dev/random_search_runoutSim_individual_performance.Rd")

## Calculate runout error ####
runout_polygons$sim_roc <- sapply(results_df, function(x) x$roc)
runout_polygons$sim_length_relerror <- sapply(results_df, function(x) x$length.relerr)
runout_polygons$sim_length_error <- sapply(results_df, function(x) x$length.error)
runout_polygons$sim_length <- runout_polygons$length + runout_polygons$sim_length_error 

median(runout_polygons$sim_length_relerror)
IQR(runout_polygons$sim_length_relerror)

median(runout_polygons$sim_length_error)


png(filename="C:/GitProjects/runoutSim/Case_Study/Figures/random_search_hist_error_plots.png", res = 300, width = 7.5, height = 6,
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

library(patchwork)  # or use gridExtra if you prefer


# Run model from training source ###############################################

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


# Define number of cores to use
n_cores <- detectCores() -2

dem_terra <- terra::rast('C:/GitProjects/runoutSim/Data/elev_fillsinks_WangLiu.tif')

packed_dem <- wrap(dem_terra)

# Create parallel loop
cl <- makeCluster(n_cores, type = "PSOCK") # Open clusters


clusterExport(cl, varlist = c("packed_dem", "best_params"))

# Load required packages to each cluster
clusterEvalQ(cl, {
  library(terra)
  library(runoutSim)
})


multi_sim_paths <- parLapply(cl, source_l, function(x) {
  
  runoutSim(dem = unwrap(packed_dem), xy = x, 
            mu = best_params$pcm_mu, 
            md = best_params$pcm_md, 
            slp_thresh = best_params$rw_slp, 
            exp_div = best_params$rw_ex, 
            per_fct = best_params$rw_per, 
            walks = 1000)
})



stopCluster(cl) 

# Coerce runout modelling to raster ############################################

trav_freq <- walksToRaster(multi_sim_paths, dem_terra)
trav_prob <- runoutSim::rasterCdf(trav_freq)
trav_vel <- velocityToRaster(multi_sim_paths, dem_terra)

## Visualize runout modelling results ####
leafmap(runout_polygons, label = "Runout observations", opacity = .35) %>%
  leafmap(trav_prob, palette =viridis::viridis(10, direction = -1), label = "Traverse probability") %>%
  leafmap(trav_vel, palette = 'plasma', label = "Velocity") %>%
  leafmap(source_points, color = "red", label = "Runout source") 

## Export runout modelling results (rasters) ####
writeRaster(trav_prob, filename = "C:\\GitProjects\\runoutSim\\Case_Study\\Figures\\randsearch_opt_runout_trav_ecdf_prob.tif")
writeRaster(trav_vel, filename = "C:\\GitProjects\\runoutSim\\Case_Study\\Figures\\randsearch_runout_trav_vel.tif")



