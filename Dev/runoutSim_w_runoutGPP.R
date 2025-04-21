# Development (Dev) environment for optimizing random walk runout simulations

# source("Dev/R/runoptGPP_Dev/pcm_gridsearch.R")
# source("Dev/R/runoptGPP_Dev/pcm_performance.R")
# source("Dev/R/runoptGPP_Dev/pcm_spatial_cross_validation.R")
# source("Dev/R/runoptGPP_Dev/randomwalk_gridsearch.R")
# source("Dev/R/runoptGPP_Dev/randomwalk_performance.R")
# source("Dev/R/runoptGPP_Dev/randomwalk_spatial_cross_validation.R")
# source("Dev/R/runoptGPP_Dev/raster_rescale.R")
# source("Dev/R/runoptGPP_Dev/runout_geometry.R")
# source("Dev/R/runoptGPP_Dev/source_area_threshold.R")

library(runoptGPP)
library(runoutSim)

# Optimizing an individual runout path simulation ##############################

library(raster)
library(sf)

# Load digital elevation model (DEM)
dem <- raster("Dev/Data/elev_nosinks.tif") # use sink filled DEM to remove pits and flats 

# Load runout source points and polygons
source_points <- st_read("Dev/Data/debris_flow_source_points.shp")
runout_polygons <- st_make_valid(st_read("Dev/Data/debris_flow_runout_polygons.shp"))
# ^ Need to clean this up so make valid not needed

# Select a single debris flow and source point for the example
runout_polygon <- runout_polygons[10,]

# Get corresponding source point
source_point  <- st_filter(st_as_sf(source_points), st_as_sf(runout_polygon))

plot(as_Spatial(runout_polygon))
plot(as_Spatial(source_point), add = TRUE)


## Random walk simulation (path) ####

rwPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
              slp = 30, ex = 3, per = 2,
              gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
              plot_eval = TRUE)

steps <- 3
# Exponent controlling lateral spread
rwexp_vec <- seq(1.3, 3, len=steps)
# Persistence factor to weight flow direction consistency
rwper_vec <- seq(1.5, 2, len=steps)
# Slope threshold - below lateral spreading is modeled.
rwslp_vec <- seq(20, 40, len=steps)

rw_gridsearch <- rwGridsearch(dem, slide_plys = runout_polygon, slide_src = source_point,
                              #Input random walk grid search space
                              slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
                              #Set number of simulation iterations
                              gpp_iter = 1000,
                              #Define processing extent size (m)
                              buffer_ext = 0,
                              #(Optional) Define size of buffer to make source area from point
                              buffer_source = 50,
                              saga_lib = NULL)

rw_gridsearch
rw_opt_single <- rwGetOpt_single(rw_gridsearch)
rw_opt_single

## PCM runout (length) ####
pcm <- pcmPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
                      rw_slp = 40, rw_ex = 2.15, rw_per = 1.5,
                      pcm_mu = 0.08, pcm_md = 40,
                      gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                      plot_eval = TRUE, return_features = TRUE)

pcm$length.error

# The mass-to-drag ratio (m)
pcmmd_vec <- seq(20, 120, by=20)
# The sliding friction coefficient
pcmmu_vec <- seq(0.05, 0.3, by=0.1)

pcm_gridsearch <- pcmGridsearch(dem,
                                slide_plys = runout_polygon, slide_src = source_point,
                                #Plug-in random walk optimal parameters
                                rw_slp = rw_opt_single$rw_slp_opt,
                                rw_ex = rw_opt_single$rw_exp_opt,
                                rw_per = rw_opt_single$rw_per_opt,
                                #Input PCM grid search space
                                pcm_mu_v = pcmmu_vec,
                                pcm_md_v = pcmmd_vec,
                                #Set number of simulation iterations
                                gpp_iter = 1000,
                                #Define processing extent size (m)
                                buffer_ext = 500,
                                #(Optional) Define size of buffer to make source area from point
                                buffer_source = 50, 
                                saga_lib = NULL)

# Get optimal parameters
pcmGetOpt_single(pcm_gridsearch)

# Regionally optimizing and spatially validating ###############################

plot(terra::rast(dem))
plot(as_Spatial(runout_polygons), add = TRUE)

steps <- 11
rwexp_vec <- seq(1.3, 3, len=steps)
rwper_vec <- seq(1.5, 2, len=steps)
rwslp_vec <- seq(20, 40, len=steps)

rwexp_vec


library(foreach)

# Define which runout polygons are used for optimization
polyid_vec <- 1:10

# Set up cluster
cl <- parallel::makeCluster(7)
doParallel::registerDoParallel(cl)

# Run grid search loop
rw_gridsearch_multi <-
  foreach(poly_id=polyid_vec, .packages=c('raster', 'ROCR', 'Rsagacmd', 'sf', 'runoutSim', 'terra')) %dopar% {
    
    #.GlobalEnv$saga <- NULL
    
    rwGridsearch(dem, slide_plys = runout_polygons, slide_src = source_points,
                 slide_id = poly_id, slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
                 gpp_iter = 1000, buffer_ext = 500, buffer_source = 50, save_res = FALSE,
                 plot_eval = FALSE, saga_lib = NULL)
    
  }

parallel::stopCluster(cl)
