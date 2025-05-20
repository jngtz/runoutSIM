# Development (Dev) environment for optimizing random walk runout simulations
library(runoutSim)

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
# source("./Dev/R/pcm.R")
# source("./Dev/R/random_walk.R")
# source("./Dev/R/simulation_to_raster.R")
# source("./Dev/R/runout_connectivity.R")
# source("./Dev/R/interactive_plot.R")

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
dem <- raster("Dev/Data/elev_nosinks.tif") # use sink filled DEM to remove pits and flats 

# Load runout source points and polygons
source_points <- st_read("Dev/Data/debris_flow_source_points.shp")
runout_polygons <- st_make_valid(st_read("Dev/Data/debris_flow_runout_polygons.shp"))
# ^ Need to clean this up so make valid not needed

# Select a single debris flow and source point for the example
#runout_polygons <- runout_polygons[10:20,]
#runout_polygons$sim_id <- 1:nrow(runout_polygons)

rungeom <- runoptGPP::runoutGeom(as_Spatial(runout_polygons), dem)
rungeom$fid <- rungeom$id <- NULL

runout_polygons <- cbind(runout_polygons, rungeom)

# Get corresponding source point
source_points  <- st_filter(st_as_sf(source_points), st_as_sf(runout_polygons))

leafmap(runout_polygons) %>% leafmap(source_points, col = "red")

# 

#library(mlrMBO) # 
#library(mlr)
#library(DiceKriging) # used for our proxy/surrogate model for BO

#library(lhs)


## Set-up working environment to allow for parallel processing #################
library(doParallel)
library(foreach)
library(ParamHelpers)
library(lhs)
library(mlr)
library(DiceKriging)
library(mlrMBO)

# Define number of cores 
ncores <- parallel::detectCores() - 3
cl     <- makeCluster(ncores)
registerDoParallel(cl)

# Make sure each worker (core) has access to the data and functions
clusterExport(cl, c("dem", "runout_polygons", "source_points", "pcmPerformance",
                    "runoutSim"))

# Load required packages on each worker
clusterEvalQ(cl, {
  library(raster); library(terra)
  library(sp);     library(sf)
  library(ROCR);   library(Rsagacmd)
  library(runoutSim); library(runoptGPP)
  TRUE
})

## Set-up parameter space, surrogate model and acquisition function ############


# Define parameters' ranges to explore in optimization. These have been set to
# approximate parameter ranges by Wichmann (2017 in NHESS, Table 3).

ps <- makeParamSet(
  makeNumericParam("slp", lower = 20, upper = 40),
  makeNumericParam("ex",  lower = 1.3,  upper = 3),
  makeNumericParam("per", lower = 1.5, upper = 2),
  makeNumericParam("mu",  lower = 0.04, upper = 0.8),
  makeNumericParam("md",  lower = 20, upper = 150)
)

# Create a learner object / define the surrogate model for optimization. In this
# case we are using kriging regression from the DiceKriging package.


lrn  <- makeLearner("regr.km", predict.type = "se")

# Create a control object to configure the behaviour of the Bayesian optimization
# process


ctrl <- makeMBOControl()

# Set the termination criteria (e.g. run for 30 iterations)
ctrl <- setMBOControlTermination(ctrl, iters = 30)

# Set the acquisition function (infill criteria). This determines where to sample
# next (in parameter space) based on the surrogate model predictions.
# makeBBOInfillCritEI() picks points that maximize the expected improvement (EI) 
# and balances exploration (trying new spaces) and exploitation (refining the 
# known good areas)

ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())

## Define a global (multi) objective function ##################################

# Here we are creating an objective function that finds optimal parameters for
# regional runout modelling. I.e. the best set of optimal parameters to reduce
# runout simulation error for the entire runout sample.

#library(smoof)

obj.fun <- makeSingleObjectiveFunction(
  name = "Global Runout Optimization",
  fn   = function(x) {
    x <- as.list(x)
    
    # We use an inner for loop (parallelized) to evaluate performance on 
    # every slide when sampling for optimization function
    
    results_df <- foreach(
      i = seq_len(nrow(runout_polygons)),
      .combine = rbind
    ) %dopar% {
      
      # runoptGPP::pcmPerformance function is used to determine error:
      # area under the roc, and runout length (bounding box method)
      # for each runout sample simulation
      
      out <- tryCatch({
        
        pcmPerformance(
          dem            = dem,
          slide_plys     = runout_polygons,
          slide_src      = source_points,
          slide_id       = i,
          rw_slp         = x$slp, # x will be pulled from our ps object
          rw_ex          = x$ex,
          rw_per         = x$per,
          pcm_mu         = x$mu,
          pcm_md         = x$md,
          gpp_iter       = 1000,
          buffer_ext     = NULL,
          buffer_source  = 20,
          plot_eval      = FALSE,
          return_features= FALSE
        )
        
      }, error = function(e) {
        message("Slide ", i, " failed: ", e$message)
        return(NULL)
      })
      
      if (is.null(out)){
        
        return(data.frame(
          id = i, 
          roc = NA, 
          len_err = NA, 
          len_rerr = NA)) # could treat these a no run possible? e.g. roc  = 0, len_err = length of polygon, len_rerr = 1
        
      } else {
        
        return(data.frame(
          id      = i, # row id of runout sample
          roc     = out$roc,
          len_err = out$length.error,
          len_rerr = out$ length.relerr,
          stringsAsFactors = FALSE))
        
      }
      
    }
    
    # We can improve this by having a training and testing error (e.g. cross 
    # validation)
    
    # remove NAs
    valid <- complete.cases(results_df)
    results_df <- results_df[valid,]
    
    
    # compute the global loss: mean(–ROC) + mean(|len_err|)
    # Need to score leng_error - if leng_err need to make from 0 to 1
    score <- 0.1*(1-mean(results_df$roc)) + 0.9*mean(results_df$len_rerr)
  },
  par.set  = ps,
  minimize = TRUE
)

## Initialize and Run Bayesian Optimization ####################################

design <- generateDesign(
  n       = 5 * length(ps$pars), # number of initial samples (5 x number of parameters)
  par.set = ps, # parameter set
  fun     = lhs::maximinLHS # latin hypercube sample (stratified sampling in parameter space)
)

# Run optimization

start_time <- Sys.time()

global_run <- mbo(
  fun      = obj.fun,
  design   = design,
  learner  = lrn,
  control  = ctrl,
  show.info= TRUE
)

end_time <- Sys.time()
run_time <- end_time - start_time
run_time #2.25 hours

## Summarize parameter optimization ############################################

cat("Best global parameters:\n")
print(global_run$x)
cat("Global objective value:\n")
print(global_run$y)


# Save results
setwd("C:\\sda\\Workspace\\sedconnect")
save(global_run, file = "Global_MBO_RunoutSim.Rd")
(load("Global_MBO_RunoutSim.Rd"))

# Check optimization if it was sufficient
plot(global_run)

# Retrieve overall performance

results_df <- foreach(
  i = seq_len(nrow(runout_polygons))) %dopar% {
    out <- 
      pcmPerformance(
        dem            = dem,
        slide_plys     = runout_polygons,
        slide_src      = source_points,
        slide_id       = i,
        rw_slp         = global_run$x$slp,
        rw_ex          = global_run$x$ex,
        rw_per         = global_run$x$per,
        pcm_mu         = global_run$x$mu,
        pcm_md         = global_run$x$md,,
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

dem_terra <- terra::rast('C:\\sda\\GitProjects\\runoutSim\\Dev\\Data\\elev_nosinks.tif')

packed_dem <- wrap(dem_terra)

# Create parallel loop
cl <- makeCluster(n_cores, type = "PSOCK") # Open clusters

# Load objects and 'custom' functions to each cluster
# clusterExport(cl, varlist = c("runoutSim", "euclideanDistance", "adjCells",
#                              "adjRowCol", "pcm", "packed_dem","sourceConnect",
#                              "feature_mask"))

clusterExport(cl, varlist = c("packed_dem", "global_run"))

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
            mu = global_run$x$mu, 
            md = global_run$x$md, 
            slp_thresh = global_run$x$slp, 
            exp_div = global_run$x$ex, 
            per_fct = global_run$x$per, 
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