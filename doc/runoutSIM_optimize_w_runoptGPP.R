## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")
library(here)
knitr::opts_knit$set(root.dir = here())

## ----eval=FALSE---------------------------------------------------------------
# remotes::install_github("jngtz/runoptGPP")

## ----message=FALSE, results='hide', fig.show='hold', fig.height=7, fig.width = 5, fig.align='center'----
# load packages
library(runoutSim)
library(runoptGPP)
library(terra)
library(sf)

# Load digital elevation model (DEM)
dem <- rast("Data/elev_fillsinks_WangLiu.tif")

# Compute hillshade for visualization 
slope <- terrain(dem, "slope", unit="radians")
aspect <- terrain(dem, "aspect", unit="radians")
hill <- round(shade(slope, aspect, 40, 270, normalize = TRUE))

# Load debris flow runout source points and polygons
source_points <- st_read("Data/debris_flow_source_points.shp")
runout_polygons <- st_read("Data/debris_flow_runout_polygons.shp")

# Plot input data
plot(hill, col=grey(150:255/255), legend=FALSE,
     mar=c(2,2,1,4))
plot(dem, col=viridis::mako(100), alpha = .5, add = TRUE)
plot(st_geometry(runout_polygons), add = TRUE)

## -----------------------------------------------------------------------------
steps <- 5
rwexp_vec <- seq(1.3, 3, len=steps) # Expondent of divergence
rwper_vec <- seq(1.5, 2, len=steps) # Persistence factor
rwslp_vec <- seq(20, 40, len=steps) # Slope threshold

rwexp_vec
rwper_vec
rwslp_vec

## ----fig.show='hold', fig.height=4, fig.width = 7, fig.align='center', eval = FALSE----
# library(foreach)
# library(raster)
# polyid_vec <- 1:nrow(source_points)
# 
# n_cores <- parallel::detectCores() -2
# cl <- parallel::makeCluster(n_cores)
# doParallel::registerDoParallel(cl)
# 
# #Coerce dem to raster() dem.
# dem <- raster(dem)
# 
# rw_gridsearch_multi  <-
#   foreach(poly_id=polyid_vec, .packages=c('terra','raster', 'ROCR', 'sf', 'runoptGPP', 'runoutSim')) %dopar% {
# 
#     rwGridsearch(dem, slide_plys = runout_polygons, slide_src = source_points,
#                  slide_id = poly_id, slp_v = rwslp_vec, ex_v = rwexp_vec,
#                  per_v = rwper_vec, gpp_iter = 1000, buffer_ext = 500, buffer_source = NULL,
#                  save_res = TRUE, plot_eval = FALSE)
# 
#   }
# 
# parallel::stopCluster(cl)

## ----include = FALSE----------------------------------------------------------
load("Data\\rw_gridsearch_multi.Rd")

## ----warning=FALSE------------------------------------------------------------
rw_opt <- rwGetOpt(rw_gridsearch_multi, 
                   measure = median)
rw_opt

## -----------------------------------------------------------------------------
# Define PCM model grid seach space
pcmmd_vec <- seq(20, 140, by=20) # mass-to-drag ratio (m)
pcmmu_vec <- seq(0.01, 0.1, by=0.01) # sliding friction coefficient 

pcmmd_vec

## ----eval = FALSE-------------------------------------------------------------
# # Run using parallelization
# cl <- parallel::detectCores() -2
# doParallel::registerDoParallel(cl)
# 
# pcm_gridsearch_multi <-
#   foreach(poly_id=polyid_vec, .packages=c('terra','raster', 'ROCR', 'sf', 'runoptGPP', 'runoutSim')) %dopar% {
# 
#     pcmGridsearch(dem,
#                   slide_plys = runout_polygons, slide_src = source_points,
#                   slide_id = poly_id, rw_slp = 40, rw_ex = 3, rw_per = 1.9,
#                   pcm_mu_v = pcmmu_vec, pcm_md_v = pcmmd_vec,
#                   gpp_iter = 1000,
#                   buffer_ext = NULL, buffer_source = NULL,
#                   predict_threshold = 0.5, save_res = FALSE)
#   }
# 
# parallel::stopCluster(cl)

## ----include = FALSE----------------------------------------------------------
load("Data\\pcm_gridsearch_multi.Rd")

## ----warning=FALSE, fig.show='hold', fig.height=4, fig.width = 7, fig.align='center'----
pcmGetOpt(pcm_gridsearch_multi, 
          performance = "relerr", 
          measure = "median", 
          plot_opt = TRUE)

## ----fig.show='hold', fig.height=4, fig.width = 7, fig.align='center', warning = FALSE----
par(mfrow = c(1,2))
rw_spcv <- rwSPCV(x = rw_gridsearch_multi, slide_plys = runout_polygons,
                  n_folds = 5, repetitions = 10)

freq_rw <- rwPoolSPCV(rw_spcv, plot_freq = TRUE)

pcm_spcv <- pcmSPCV(pcm_gridsearch_multi, slide_plys = runout_polygons,
                    n_folds = 5, repetitions = 10, from_save = FALSE)

freq_pcm <- pcmPoolSPCV(pcm_spcv, plot_freq = TRUE)

freq_rw
freq_pcm

