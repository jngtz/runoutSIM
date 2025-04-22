# Development (Dev) environment for optimizing random walk runout simulations

source("Dev/R/runoptGPP_Dev/pcm_gridsearch.R")
source("Dev/R/runoptGPP_Dev/pcm_performance.R")
source("Dev/R/runoptGPP_Dev/pcm_spatial_cross_validation.R")
source("Dev/R/runoptGPP_Dev/randomwalk_gridsearch.R")
source("Dev/R/runoptGPP_Dev/randomwalk_performance.R")
source("Dev/R/runoptGPP_Dev/randomwalk_spatial_cross_validation.R")
source("Dev/R/runoptGPP_Dev/raster_rescale.R")
source("Dev/R/runoptGPP_Dev/runout_geometry.R")
source("Dev/R/runoptGPP_Dev/source_area_threshold.R")

#library(runoutSim)
source("./Dev/R/pcm.R")
source("./Dev/R/random_walk.R")
source("./Dev/R/simulation_to_raster.R")
source("./Dev/R/runout_connectivity.R")
source("./Dev/R/interactive_plot.R")

# New strategy
# - 1st estimate runout distance PCM (quicker)
# - 2nd estimate runout path (shorter runout will help reduce processing times.)

# Task 04.22.2025
# - transfer update in pcm_performance.R (runoptGPP)

# Optimizing an individual runout path simulation ##############################

library(raster)
library(terra)
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
pcm_runout <- pcmPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
                             rw_slp = 40, rw_ex = 2.15, rw_per = 1.5,
                             pcm_mu = 0.08, pcm_md = 40,
                             gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                             plot_eval = TRUE, return_features = TRUE)

pcm_runout$length.error

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

setwd("C:\\sda\\Workspace\\sedconnect")


plot(terra::rast(dem))
plot(as_Spatial(runout_polygons), add = TRUE)

steps <- 11
rwexp_vec <- seq(1.3, 3, len=steps)
rwper_vec <- seq(1.5, 2, len=steps)
rwslp_vec <- seq(20, 40, len=steps)

rwexp_vec


library(foreach)

# Define which runout polygons are used for optimization
polyid_vec <- 1:nrow(runout_polygons)

# Set up cluster
cl <- parallel::makeCluster(28)
doParallel::registerDoParallel(cl)

# Run grid search loop
rw_gridsearch_multi <-
  foreach(poly_id=polyid_vec, .packages=c('raster', 'ROCR', 'Rsagacmd', 'sf', 'runoutSim', 'terra')) %dopar% {
    
    #.GlobalEnv$saga <- NULL
    
    rwGridsearch(dem, slide_plys = runout_polygons, slide_src = source_points,
                 slide_id = poly_id, slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
                 gpp_iter = 1000, buffer_ext = 500, buffer_source = 50, save_res = TRUE,
                 plot_eval = FALSE, saga_lib = NULL)
    
  }

parallel::stopCluster(cl)


pcmmd_vec <- seq(20, 150, by=5) # mass-to-drag ratio (m)
pcmmu_vec <- seq(0.04, 0.6, by=0.01) # sliding friction coefficient 

pcmmd_vec

# Run using parallelization
cl <- parallel::makeCluster(28)
doParallel::registerDoParallel(cl)

pcm_gridsearch_multi <-
  foreach(poly_id=polyid_vec, .packages=c('raster', 'ROCR', 'Rsagacmd', 'sf', 'runoutSim', 'terra')) %dopar% {
    
    .GlobalEnv$saga <- saga
    
    pcmGridsearch(dem,
                  slide_plys = runout_polygons, slide_src = source_points, slide_id = poly_id,
                  rw_slp = rw_opt$rw_slp_opt, rw_ex = rw_opt$rw_exp_opt, rw_per = rw_opt$rw_per_opt,
                  pcm_mu_v = pcmmu_vec, pcm_md_v = pcmmd_vec,
                  gpp_iter = 1000,
                  buffer_ext = 500, buffer_source = NULL,
                  predict_threshold = 0.5,
                  plot_eval = FALSE, saga_lib = saga)
    
  }

parallel::stopCluster(cl)


# Test multi-objective Baysiean optimization ###################################

library(mlrMBO)
library(mlr)
library(DiceKriging)
library(ParamHelpers)

# Define the parameter space
ps <- makeParamSet(
  makeIntegerParam("slp", lower = 20L, upper = 50L),   # slope threshold
  makeIntegerParam("ex", lower = 1L, upper = 3L),      # lateral spread exponent
  makeNumericParam("per", lower = 1.5, upper = 2),     # persistence
  makeNumericParam("mu", lower = 0.04, upper = 0.6)    # friction
  # md is fixed in the function
)

slide <- 10

# Updated objective function
obj.fun <- makeSingleObjectiveFunction(
  name = "My Runout Optimization",
  fn = function(x) {
    x <- as.list(x)    # <— now x$slp, x$ex, etc. will work
    res <- pcmPerformance(
      dem          = dem,
      slide_plys   = runout_polygons,
      slide_src    = source_points,
      slide_id     = slide,
      rw_slp       = x$slp,
      rw_ex        = x$ex,
      rw_per       = x$per,
      pcm_mu       = x$mu,
      pcm_md       = 40,
      gpp_iter     = 1000,
      buffer_ext   = 500,
      buffer_source= 20,
      plot_eval    = FALSE,
      return_features = FALSE
    )
    # Note: pcmPerformance returns `res$roc` and `res$length.error`
    - res$roc + abs(res$length.error)
  },
  par.set  = ps,
  minimize = TRUE
)

# Kriging surrogate learner from mlr
lrn <- makeLearner("regr.km", predict.type = "se")

# Configure MBO
ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters = 20)
ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())

# Run optimization
design <- generateDesign(n = 5 * length(ps$pars),
                         par.set = ps,
                         fun     = lhs::maximinLHS)

run <- mbo(obj.fun,
           design  = design,
           learner = lrn,
           control = ctrl)

plot(run)

# plot optimal
run$x

pcm_runout <- pcmPerformance(dem, slide_plys = runout_polygons, slide_src = source_points,
                             slide_id = slide,
                             rw_slp = run$x$slp, rw_ex = run$x$ex, rw_per = run$x$per,
                             pcm_mu = run$x$mu, pcm_md = 40,
                             gpp_iter = 1000, buffer_ext = 500, buffer_source = 40,
                             plot_eval = TRUE, return_features = TRUE)

run_res <- terra::rast(pcm_runout$gpp.parea)
ply_res <- st_as_sf(pcm_runout$actual.poly)
leafplot(data = run_res) %>% leafplot(data = ply_res)


# Test regional optimization ####################################################

## Test runout on all slides

for(i in 29:nrow(runout_polygons)){
  print(i)
  pcm_runout <- pcmPerformance(dem, slide_plys = runout_polygons, slide_src = source_points,
                               slide_id = i,
                               rw_slp = run$x$slp, rw_ex = run$x$ex, rw_per = run$x$per,
                               pcm_mu = run$x$mu, pcm_md = 40,
                               gpp_iter = 1000, buffer_ext = 1000, buffer_source = 40,
                               plot_eval = TRUE, return_features = TRUE)
}


library(mlrMBO)
library(mlr)
library(DiceKriging)
library(ParamHelpers)
library(lhs)
library(doParallel)
library(foreach)

# 1) Set up your fixed objects
ps <- makeParamSet(
  makeIntegerParam("slp", lower = 20L, upper = 50L),
  makeIntegerParam("ex",  lower = 1L, upper = 3L),
  makeNumericParam("per", lower = 1.5, upper = 2),
  makeNumericParam("mu",  lower = 0.04, upper = 0.6)
)

lrn <- makeLearner("regr.km", predict.type = "se")
ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters = 20)
ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())

# 2) Parallel backend
Sys.time()
n.cores <- parallel::detectCores() - 1
cl <- makeCluster(20)
registerDoParallel(cl)

# 3) Slide IDs
slide_ids <- 1:nrow(runout_polygons)

# 4) Run MBO for each slide in parallel
results <- foreach(
  slide = slide_ids,
  .packages = c("mlrMBO","mlr","DiceKriging","ParamHelpers","lhs",
                'raster', 'ROCR', 'Rsagacmd', 'sf', 'runoutSim', 'terra'),
  .combine  = rbind
) %dopar% {
  
  # build objective for this slide
  obj.fun <- makeSingleObjectiveFunction(
    name = paste0("runout_slide_", slide),
    fn = function(x) {
      x <- as.list(x)
      res <- pcmPerformance(
        dem          = dem,
        slide_plys   = runout_polygons,
        slide_src    = source_points,
        slide_id     = slide,
        rw_slp       = x$slp,
        rw_ex        = x$ex,
        rw_per       = x$per,
        pcm_mu       = x$mu,
        pcm_md       = 40,
        gpp_iter     = 1000,
        buffer_ext   = 500,
        buffer_source= 20,
        plot_eval    = FALSE,
        return_features = FALSE
      )
      # minimize: –AUC + |length.err|
      return(-res$roc + abs(res$length.error))
    },
    par.set  = ps,
    minimize = TRUE
  )
  
  # initial design
  design <- generateDesign(
    n       = 5 * length(ps$pars),
    par.set = ps,
    fun     = maximinLHS
  )
  
  # run MBO
  run <- mbo(
    fun     = obj.fun,
    design  = design,
    learner = lrn,
    control = ctrl,
    show.info = FALSE
  )
  
  # return a one‐row data.frame of results
  data.frame(
    slide = slide,
    slp   = run$x$slp,
    ex    = run$x$ex,
    per   = run$x$per,
    mu    = run$x$mu,
    obj   = run$y,
    stringsAsFactors = FALSE
  )
}

stopCluster(cl)

# results is now a data.frame with one row per slide
print(results)
Sys.time()