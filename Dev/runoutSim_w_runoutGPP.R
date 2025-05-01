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

rwexp_vec <- seq(1.2, 3, by=.2)  
rwper_vec <- seq(1.5, 2, by=0.05)   
rwslp_vec <- seq(25, 35, by=1)   
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
save(rw_gridsearch_multi, file = "rw_gridsearch.Rd")


setwd("C:\\sda\\Workspace\\sedconnect")
(load("rw_gridsearch.Rd"))
rw_opt <- rwGetOpt(rw_gridsearch_multi, measure = median)
rw_opt

pcmmd_vec <- seq(20, 150, by=5) # mass-to-drag ratio (m)
pcmmu_vec <- seq(0.04, 0.6, by=0.01) # sliding friction coefficient 

pcmmd_vec

# Run using parallelization
cl <- parallel::makeCluster(28)
doParallel::registerDoParallel(cl)

pcm_gridsearch_multi <-
  foreach(poly_id=polyid_vec, .packages=c('raster', 'ROCR', 'Rsagacmd', 'sf', 'runoutSim', 'terra')) %dopar% {
    
    #.GlobalEnv$saga <- saga
    
    pcmGridsearch(dem,
                  slide_plys = runout_polygons, slide_src = source_points, slide_id = poly_id,
                  rw_slp = rw_opt$rw_slp_opt, rw_ex = rw_opt$rw_exp_opt, rw_per = rw_opt$rw_per_opt,
                  pcm_mu_v = pcmmu_vec, pcm_md_v = pcmmd_vec,
                  gpp_iter = 1000,
                  buffer_ext = 500, buffer_source = NULL,
                  predict_threshold = 0.5,
                  save_res = TRUE,
                  plot_eval = FALSE, saga_lib = NULL)
    
  }

parallel::stopCluster(cl)


# Test multi-objective Baysiean optimization ###################################

# Would be good for hybrid modelling - using custom cost function 

library(mlrMBO)
library(mlr)
library(DiceKriging)
library(ParamHelpers)

# Define the parameter space
ps <- makeParamSet(
  makeIntegerParam("slp", lower = 20, upper = 50),   # slope threshold
  makeIntegerParam("ex", lower = 1, upper = 3),      # lateral spread exponent
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
# This optimizes each landslide individually

## Test runout on all slides

for(i in 60:nrow(runout_polygons)){
  print(i)
  pcm_runout <- pcmPerformance(dem, slide_plys = runout_polygons, slide_src = source_points,
                               slide_id = i,
                               rw_slp = run$x$slp, rw_ex = run$x$ex, rw_per = run$x$per,
                               pcm_mu = run$x$mu, pcm_md = 40,
                               gpp_iter = 1000, buffer_ext = 1000, buffer_source = 30,
                               plot_eval = TRUE, return_features = TRUE)
}

# ────────────────────────────────────────────────────────────────────────────
# 1) Libraries & cluster setup
# ────────────────────────────────────────────────────────────────────────────

# Not working, error is in the foreah loop setup, need to better put it together...

library(mlrMBO)
library(mlr)
library(DiceKriging)
library(ParamHelpers)
library(lhs)
library(doParallel)
library(foreach)

# number of cores for the inner foreach
ncores <- parallel::detectCores() - 1
cl     <- makeCluster(ncores)
registerDoParallel(cl)

# make sure each worker has access to your data + function
clusterExport(cl, c("dem", "runout_polygons", "source_points", "pcmPerformance",
                    "runoutSim"))

# load required packages on each worker
clusterEvalQ(cl, {
  library(raster); library(terra)
  library(sp);     library(sf)
  library(ROCR);   library(Rsagacmd)
  library(runoutSim); library(runoptGPP)
  TRUE
})

# ────────────────────────────────────────────────────────────────────────────
# 2) Parameter space, surrogate & MBO control
# ────────────────────────────────────────────────────────────────────────────
ps <- makeParamSet(
  makeIntegerParam("slp", lower = 20L, upper = 50L),
  makeIntegerParam("ex",  lower = 1L,  upper = 3L),
  makeNumericParam("per", lower = 1.5, upper = 2),
  makeNumericParam("mu",  lower = 0.04, upper = 0.6)
)

lrn  <- makeLearner("regr.km", predict.type = "se")
ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters = 30)
ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())

# ────────────────────────────────────────────────────────────────────────────
# 3) Global objective function
# ────────────────────────────────────────────────────────────────────────────
obj.fun <- makeSingleObjectiveFunction(
  name = "Global Runout Optimization",
  fn   = function(x) {
    x <- as.list(x)
    
    # inner parallel loop: evaluate pcmPerformance on every slide
    results_df <- foreach(
      i = seq_len(nrow(runout_polygons)),
      .combine = rbind
    ) %dopar% {
      out <- 
        pcmPerformance(
          dem            = dem,
          slide_plys     = runout_polygons,
          slide_src      = source_points,
          slide_id       = i,
          rw_slp         = x$slp,
          rw_ex          = x$ex,
          rw_per         = x$per,
          pcm_mu         = x$mu,
          pcm_md         = 40,
          gpp_iter       = 1000,
          buffer_ext     = 500,
          buffer_source  = 20,
          plot_eval      = FALSE,
          return_features= FALSE
        )

      data.frame(
        id      = i,
        roc     = out$roc,
        len_err = out$length.error,
        #len_rerr = out$ length.relerr,
        stringsAsFactors = FALSE
      )
    }
    
    # remove NAs
    valid <- complete.cases(results_df)
    results_df <- results_df[valid,]

    
    # compute the global loss: mean(–ROC) + mean(|len_err|)
    score <- -mean(results_df$roc) + mean(abs(results_df$len_err))
  },
  par.set  = ps,
  minimize = TRUE
)

# ────────────────────────────────────────────────────────────────────────────
# 4) Initial design & run MBO
# ────────────────────────────────────────────────────────────────────────────
design <- generateDesign(
  n       = 5 * length(ps$pars),
  par.set = ps,
  fun     = lhs::maximinLHS
)

global_run <- mbo(
  fun      = obj.fun,
  design   = design,
  learner  = lrn,
  control  = ctrl,
  show.info= TRUE
)

# ────────────────────────────────────────────────────────────────────────────
# 5) Clean up & results
# ────────────────────────────────────────────────────────────────────────────
stopCluster(cl)

cat("Best global parameters:\n")
print(global_run$x)
cat("Global objective value:\n")
print(global_run$y)

# PRoposed approach ###########################################################

# We use Bayesian optimization to find the optimal values of a spatial mu
# that reduce the runout length

# the random walk parameters will be estimated on experimentation or optimized
# without pcm, then we will focus on optimizing only pcm



#https://medium.com/data-science/bayesian-optimization-concept-explained-in-layman-terms-1d2bcdeaf12f


# Libraries (load or install)
library(mlrMBO)
library(ParamHelpers)
library(mlr)  # Required for some mlrMBO dependencies
library(ggplot2)

# Your raster covariates (e.g., forest, carea, precip) must be pre-loaded and aligned rasters
# Your fixed parameters for PCM (e.g., fix_slp, fix_ex, fix_md, etc.) should be defined
# Your pcmPerformance() function should be in memory

# Parameter space definition
ps <- makeParamSet(
  makeNumericParam("intercept", lower = -10, upper = 10),
  makeNumericParam("x1", lower = 2,  upper = 5),
  makeNumericParam("x2", lower = 1.5, upper = 2),
  makeNumericParam("x3", lower = 0.04, upper = 0.6)
)

# Objective function
sp_obj.fun <- makeSingleObjectiveFunction(
  name = "Spatial mu Optimization",
  fn = function(x) {
    # Extract values
    intercept <- x$intercept
    x1 <- x$x1
    x2 <- x$x2
    x3 <- x$x3
    
    # Build spatial mu from covariates
    sp_mu <- intercept + x1 * forest + x2 * carea + x3 * precip
    
    # may need to transform sp_mu into a range of meaninful values
    
    # Run PCM and get performance
    res <- pcmPerformance(
      dem            = dem,
      slide_plys     = runout_polygons,
      slide_src      = source_points,
      slide_id       = slide,
      rw_slp         = fix_slp,
      rw_ex          = fix_ex,
      rw_per         = fix_per,
      pcm_mu         = sp_mu,
      pcm_md         = fix_md,
      gpp_iter       = 1000,
      buffer_ext     = 500,
      buffer_source  = 20,
      plot_eval      = FALSE,
      return_features= FALSE
    )
    
    return(res$length.relerror)  # Minimize this
  },
  par.set = ps,
  minimize = TRUE
)

# Control settings
ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters = 30)  # Adjust iterations as needed
ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())  # Expected Improvement

# Surrogate model: using kriging
surrogate.lrn <- makeLearner("regr.km", predict.type = "se")

# Run Bayesian Optimization
design <- generateDesign(n = 10, par.set = ps, fun = lhs::maximinLHS)

run <- mbo(fun = sp_obj.fun, design = design, learner = surrogate.lrn, control = ctrl)

# View results
print(run)
plot(run)