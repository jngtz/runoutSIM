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
# 7. Coerce modelling results to rasters

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

# For data handling and visualization
library(ggplot2)
library(dplyr)


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

## Define grid search space ####
pcmmu_vec <- seq(0.04, 0.6, by=0.01)
polyid_vec <- 1:nrow(source_points)

## Perform grid search optimization with parallelization ####

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


# Get optimal (PCM) parameters #################################################

pcm_gridsearch_multi <- list()
files <- list.files(pattern = "result_pcm_gridsearch_")
for (i in 1:length(files)) {
  res_nm <- paste("result_pcm_gridsearch_", i, ".Rd", 
                  sep = "")
  res_obj_nm <- load(res_nm)
  result_pcm <- get(res_obj_nm)
  pcm_gridsearch_multi[[i]] <- result_pcm
}


pcm_opt <- pcmGetOpt(pcm_gridsearch_multi, performance = "relerr", 
                     measure = "median", plot_opt = FALSE, from_save = TRUE)
pcm_opt

save(pcm_opt, file = "pcm_opt_params.Rd")
save(pcm_gridsearch_multi, file = "pcm_gridsearch_multi.Rd")

# Validate model performance ###################################################

# Create the plot of median relative error to sliding friction coefficient 
# with IQR for uncertainty

grid_res <- data.frame(
  mu = as.numeric(rownames(pcmGetGrid(performance = "relerr", measure = "median", from_save = TRUE))),
  relerr_median = pcmGetGrid(performance = "relerr", measure = "median", from_save = TRUE)[,1],
  error_median = pcmGetGrid(performance = "error", measure = "median", from_save = TRUE)[,1],
  abs_error_median = abs(pcmGetGrid(performance = "error", measure = "median", from_save = TRUE)[,1]),
  iqr = pcmGetGrid(performance = "relerr", measure = "IQR", from_save = TRUE)[,1]
)


grid_res <- grid_res %>%
  mutate(
    ymin = relerr_median - iqr / 2,
    ymax = relerr_median + iqr / 2
  )

# Plot optimization results
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


pcm_spcv <- pcmSPCV(pcm_gridsearch_multi, slide_plys = runout_polygons,
                    n_folds = 5, repetitions = 100, from_save = FALSE)

freq_pcm <- pcmPoolSPCV(pcm_spcv, plot_freq = FALSE)
freq_pcm


# Retrieve overall performance

n_cores <- detectCores() -2
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

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
        pcm_mu         = as.numeric(pcm_opt$pcm_mu),
        pcm_md         = 40,
        gpp_iter       = 1000,
        buffer_ext     = NULL,
        buffer_source  = 20,
        plot_eval      = FALSE,
        return_features= TRUE
      )
  }

results_df

parallel::stopCluster(cl)

save(results_df, file = "runoutSim_individual_performance.Rd")

## Calculate runout error ####
runout_polygons$sim_roc <- sapply(results_df, function(x) x$roc)
runout_polygons$sim_length_relerror <- sapply(results_df, function(x) x$length.relerr)
runout_polygons$sim_length_error <- sapply(results_df, function(x) x$length.error)
runout_polygons$sim_length <- runout_polygons$length + runout_polygons$sim_length_error 

median(runout_polygons$sim_length_relerror)
IQR(runout_polygons$sim_length_relerror)

median(runout_polygons$sim_length_error)


png(filename="C:/GitProjects/runoutSim/Case_Study/Figures/hist_error_plots.png", res = 300, width = 7.5, height = 6,
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

# Individual plots
p1 <- ggplot(grid_res, aes(x = mu, y = relerr_median)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "gray80", alpha = 0.5) +
  geom_line(color = "#19191a", linewidth = 0.7) +
  geom_point(color = "#19191a", size = 0.7) +
  scale_x_continuous(breaks = seq(0, 2, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(
    x = expression("Sliding friction coefficient (" * mu * ")"),
    y = "Median relative runout distance error"
  ) +
  theme_light() +
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.5, "cm"),
        text = element_text(size = 7), 
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

p2 <- ggplot(runout_polygons, aes(x = sim_length_relerror)) +
  geom_histogram(bins = 40, fill = "gray70", color = "black") +
  labs(x = "Runout length relative error", y = "Frequency") +
  theme_light() +
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.5, "cm"),
        text = element_text(size = 7), 
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

p3 <- ggplot(runout_polygons, aes(x = sim_length_error)) +
  geom_histogram(bins = 40, fill = "gray70", color = "black") +
  labs(x = "Runout length error (m)", y = "Frequency") +
  theme_light() +
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.5, "cm"),
        text = element_text(size = 7), 
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

p4 <- ggplot(runout_polygons, aes(x = sim_length, y = length)) +
  geom_point(shape = 20) +
  labs(x = "Simulated runout length (m)",
       y = "Observed runout length (m)") +
  theme_light() +
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.5, "cm"),
        text = element_text(size = 7), 
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Save combined figure
perf_opt_plot <- (p1 | p2) / (p3 | p4)
perf_opt_plot <- perf_opt_plot + plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)'))) &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 9, hjust = 0, vjust = 0, face = 'bold'))
perf_opt_plot

ggsave("C:/GitProjects/runoutSim/Case_Study/Figures/hist_error_plots.png",
       plot = perf_opt_plot,
       width = 170, height = 140, units = "mm", dpi = 300)

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

dem_terra <- terra::rast('C:/GitProjects/runoutSim/Dev/Data/elev_fillsinks_WangLiu.tif')

packed_dem <- wrap(dem_terra)

# Create parallel loop
cl <- makeCluster(n_cores, type = "PSOCK") # Open clusters


clusterExport(cl, varlist = c("packed_dem"))

# Load required packages to each cluster
clusterEvalQ(cl, {
  library(terra)
  library(runoutSim)
})


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
writeRaster(trav_prob, filename = "C:\\GitProjects\\runoutSim\\Case_Study\\Results\\pcm_mu_opt_runout_trav_ecdf_prob.tif")
writeRaster(trav_vel, filename = "C:\\GitProjects\\runoutSim\\Case_Study\\Results\\pcm_mu_opt_runout_trav_vel.tif")



