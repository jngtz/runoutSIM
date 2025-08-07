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

library(mlrMBO)
library(mlr)
library(DiceKriging)
library(ParamHelpers)
library(lhs)
library(doParallel)
library(foreach)

Sys.time()
# number of cores for the inner foreach
ncores <- parallel::detectCores() - 3
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
  makeNumericParam("slp", lower = 20, upper = 50),
  makeNumericParam("ex",  lower = 1,  upper = 3),
  makeNumericParam("per", lower = 1.5, upper = 2),
  makeNumericParam("mu",  lower = 0.01, upper = 0.6),
  makeNumericParam("md",  lower = 20, upper = 50)
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
      
      out <- tryCatch({
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
      }, error = function(e) {
        message("Slide ", i, " failed: ", e$message)
        return(NULL)
      })
      
      if (is.null(out)) return(data.frame(
        id = i, roc = NA, len_err = NA, len_rerr = NA # could treat these a no run possible? e.g. roc  = 0, len_err = length of polygon, len_rerr = 1
      )) else {
        return(data.frame(
          id      = i,
          roc     = out$roc,
          len_err = out$length.error,
          len_rerr = out$ length.relerr,
          stringsAsFactors = FALSE
        ))
      }
      
      
    }
    
    # We can improve this by having a training and testing error (e.g. cross validation)
    
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

Sys.time()
# ────────────────────────────────────────────────────────────────────────────
# 5) Clean up & results
# ────────────────────────────────────────────────────────────────────────────
#stopCluster(cl)

cat("Best global parameters:\n")
print(global_run$x)
cat("Global objective value:\n")
print(global_run$y)

setwd("C:\\sda\\Workspace\\sedconnect")
save(global_run, file = "Global_MBO.Rd")
# See results scores



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
      pcm_md         = 40,
      gpp_iter       = 1000,
      buffer_ext     = 10000,
      buffer_source  = 20,
      plot_eval      = FALSE,
      return_features= TRUE
    )
}

results_df

stopCluster(cl)



# Validate

# Extract error
runout_polygons$sim_roc <- sapply(results_df, function(x) x$roc)
runout_polygons$sim_length_relerror <- sapply(results_df, function(x) x$length.relerr)
runout_polygons$sim_length_error <- sapply(results_df, function(x) x$length.error)
runout_polygons$sim_length <- runout_polygons$length + runout_polygons$sim_length_error 

plot(runout_polygons$sim_length, runout_polygons$length,
     xlab = "Simulated runout length (m)",
     ylab = "Observed runout length (m)", pch = 20)


# > change below... 


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