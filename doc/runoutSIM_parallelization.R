## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.show='hold', fig.height=7, fig.width = 5, fig.align='center'---------
# load packages
library(runoutSim)
library(terra)
library(sf)

# Load digital elevation model (DEM)
dem <- rast("C:/GitProjects/runoutSIM/Dev/Data/elev_fillsinks_WangLiu.tif")

# Compute hillshade for visualization 
slope <- terrain(dem, "slope", unit="radians")
aspect <- terrain(dem, "aspect", unit="radians")
hill <- round(shade(slope, aspect, 40, 270, normalize = TRUE))

# Load debris flow runout source points and polygons
source_points <- st_read("C:/GitProjects/runoutSIM/Dev/Data/debris_flow_source_points.shp")
runout_polygons <- st_read("C:/GitProjects/runoutSIM/Dev/Data/debris_flow_runout_polygons.shp")

# Plot input data
plot(hill, col=grey(150:255/255), legend=FALSE,
     mar=c(2,2,1,4))
plot(dem, col=viridis::mako(10), alpha = .5, add = TRUE)
plot(st_geometry(runout_polygons), add = TRUE)

# Get coordinates of source points to create a source list object
source_l <- makeSourceList(source_points)




## ----message=FALSE, results='hide'--------------------------------------------
library(parallel)

# Define number of cores to use
n_cores <- detectCores() -2

# Pack the DEM so it can be passed over a serialized connection
packed_dem <- wrap(dem)

# Create parallel loop
cl <- makeCluster(n_cores, type = "PSOCK") # Open clusters

# Export the packed DEM to each node
clusterExport(cl, varlist = c("packed_dem"))

# Load required packages to each cluster
clusterEvalQ(cl, {
  library(terra)
  library(runoutSim)
})

# Compute multiple runout simulations from source cells
multi_sim_runs <- parLapply(cl, source_l, function(x) {
  
  runoutSim(dem = unwrap(packed_dem), xy = x, 
            mu = 0.08, 
            md = 40, 
            slp_thresh = 40, 
            exp_div = 3, 
            per_fct = 1.9, 
            walks = 1000)
})
 
stopCluster(cl) 


## ----fig.show='hold', fig.height=4, fig.width = 7, fig.align='center'---------
# Coerce results to a raster
trav_freq <- walksToRaster(multi_sim_runs, method = "freq", dem)
vel_ms <- velocityToRaster(multi_sim_runs, dem, method = "max")
trav_prob <- walksToRaster(multi_sim_runs, method = "max_cdf_prob", dem)


# Plot random walks from mulitple source cells
par(mfrow = c(1,3))
plot(hill, col=grey(150:255/255), legend=FALSE, mar=c(2,2,1,4), 
     main = "Traverse frequency")
plot(trav_freq, add = TRUE, alpha = 0.7)

plot(hill, col=grey(150:255/255), legend=FALSE, mar=c(2,2,1,4), 
     main = "Traverse probability (max.)")
plot(trav_prob, add = TRUE, alpha = 0.7)

plot(hill, col=grey(150:255/255), legend=FALSE, mar=c(2,2,1,4), 
     main = "Velocity (m/s)")
plot(vel_ms, col = map.pal("plasma"), add = TRUE, alpha = 0.7)

