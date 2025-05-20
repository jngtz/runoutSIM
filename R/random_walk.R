#' Get Adjacent Cell Coordinates
#'
#' Calculates the coordinates of 8 adjacent cells surrounding a given center cell.
#'
#' @param r A numeric vector of length 2 representing the resolution of the raster in x and y directions (e.g., \code{c(xres, yres)}).
#' @param xy A numeric vector of length 2 specifying the x and y coordinates of the center cell (e.g., \code{c(x, y)}).
#'
#' @return A matrix of 8 rows and 2 columns containing the coordinates of the adjacent cells.
#' @keywords internal
#' @examples 
#' adjCells(c(10, 10), c(100, 100))
#' @export

adjCells <- function(r, xy){
  #r: vector - resolution of raster c(_,_)
  #xy: vector - xy location of center cell c(_,_)
  
  d <- c(rep(xy[1]-r[1], 3), rep(xy[1]+r[1],3), xy[1], xy[1],
         rep(c(xy[2]+r[2], xy[2], xy[2]-r[2]), 2), xy[2]+r[2], xy[2]-r[2])
  
  d <- matrix(d, ncol=2)

  
}

#' Get Row/Column Indices of Adjacent Cells
#'
#' Computes the row and column indices of the 8 adjacent cells surrounding a given central cell in a raster grid.
#'
#' @param rowcol A numeric vector of length 2 representing the row and column index of the center cell.
#'
#' @return A matrix of 8 rows and 2 columns with the row and column indices of the adjacent cells.
#' @keywords internal
#' @examples
#' adjRowCol(c(10, 10))
#' @export

adjRowCol <- function(rowcol){
  
  x <- .subset2(rowcol, 1)
  y <- .subset2(rowcol, 2)
  
  d <- c(x - 1, y -1,
         x, y -1,
         x + 1, y -1,
         
         x - 1, y + 1,
         x, y + 1,
         x + 1, y  + 1,
         
         x - 1, y,
         x + 1, y)
  
  d <- matrix(d, ncol = 2, byrow = TRUE)
  
  return(d)
}


#' Calculate Euclidean Distance Between Two Points
#'
#' Computes the Euclidean distance between two 2D points.
#'
#' @param p1 A numeric vector of length 2 representing the first point's coordinates (x, y).
#' @param p2 A numeric vector of length 2 representing the second point's coordinates (x, y).
#'
#' @return A numeric value representing the distance between the two points.
#' @examples
#' euclideanDistance(c(0, 0), c(3, 4))
#' @export

euclideanDistance <- function(p1, p2){
  sqrt( (p1[1] - p2[1])^2 + (p1[2]- p2[2])^2 )
}


#' Simulate Runout Paths Using Random Walk and PCM Physics
#'
#' Simulates the runout paths of mass movements (e.g., debris flows, snow avalanches) from one or more source points 
#' using a physically-informed random walk model. Each path incorporates slope-based directionality, 
#' lateral dispersion, and energy loss based on the Perla-Cheng-McClung (PCM) velocity model.
#'
#' @param dem SpatRaster. A digital elevation model (DEM), ideally sink-filled and hydrologically correct.
#' @param xy Numeric vector or matrix. Coordinates of one or more source points (e.g., `cbind(x, y)`).
#' @param mu Numeric. Basal friction coefficient used in the PCM velocity model (default = `0.1`).
#' @param md Numeric. Mass-to-drag ratio used in the PCM model to control deceleration from air drag and terrain resistance.
#' @param int_vel Numeric. Initial velocity for PCM model.
#' @param slp_thresh Numeric. Minimum slope (in degrees) required for particle movement. Particles stop when local slope falls below this value.
#' @param exp_div Numeric. Lateral dispersion exponent controlling spread of paths away from steepest descent. Higher values = less lateral dispersion.
#' @param per_fct Numeric. Persistence factor controlling how strongly particles follow the current downslope direction. Values >1 increase directional inertia.
#' @param walks Integer. Number of random walk particles to simulate per source point (default = `1000`).
#' @param source_connect Logical. If `TRUE`, particles are flagged based on whether they intersect a defined connectivity feature (e.g., stream channel).
#' @param connect_feature SpatRaster. Optional raster of connectivity features created using `makeConnFeature()`. Required if `source_connect = TRUE`.
#'
#' @return A list of simulated particle paths. Each element contains a matrix or list of stepwise information
#'         (e.g., coordinates, step index, velocity, slope, and connectivity).
#'
#' @details
#' This function simulates gravity-driven runout from user-defined source cells. Particle paths evolve through
#' a random walk process biased by terrain slope and lateral spread rules. At each step, the PCM velocity model
#' updates particle speed based on slope, friction (`mu`), and mass-to-drag ratio (`md`).
#'
#' Paths terminate when slope falls below `slp_thresh`, or particles encounter terrain conditions that no longer allow
#' continued downhill movement.
#'
#' @section PCM Model:
#' The PCM model (Perla-Cheng-McClung) describes velocity as a balance between gravitational acceleration,
#' frictional resistance, and air drag. Higher `md` values imply greater momentum and less deceleration from drag.
#'
#' @section Applicability:
#' This model can be applied to a variety of gravity-driven mass movements such as **debris flows**, **snow avalanches**, 
#' and other similar phenomena. The random walk approach, combined with the PCM velocity model, makes it highly versatile for 
#' simulating different types of runout based on the underlying physical principles of mass movement.
#'
#' @section Connectivity:
#' If `source_connect = TRUE`, simulated paths are checked for intersection with a connectivity feature raster.
#' This enables downstream analysis of connectivity probability (`connToRaster()`), or flow velocity surfaces (`velocityToRaster()`).
#'
#' @seealso \code{\link{pcm}}, \code{\link{makeConnFeature}}, \code{\link{walksToRaster}}, 
#'          \code{\link{connToRaster}}, \code{\link{velocityToRaster}}
#'
#' @examples
#' \dontrun{
#' Load DEM and connectivity feature
#' dem <- terra::rast("Data/elev_nosinks.tif")
#' river <- sf::st_read("Data/river_channel.shp")
#' feature_mask <- makeConnFeature(river, dem)
#'
#' # Simulate runout from a single source point
#' sim_paths <- runoutSim(dem = dem, xy = c(500000, 5200000), mu = 0.08, md = 140,
#'                        slp_thresh = 35, exp_div = 3, per_fct = 1.95, walks = 1000,
#'                        source_connect = TRUE, feature_layer = feature_mask)
#'
#' # Convert results to raster
#' trav_freq <- walksToRaster(sim_paths, dem)
#' trav_prob <- rasterCdf(trav_freq)
#' conn_prob <- connToRaster(sim_paths, dem)
#' vel_map   <- velocityToRaster(sim_paths, dem)
#' }
#' @export

runoutSim <- function(dem, xy, mu = 0.1, md = 40, int_vel = 1, slp_thresh = 30, exp_div = 3, per_fct = 2, walks = 100,
                      source_connect = FALSE, connect_feature = NULL){
  
  #---Parameter Description---
  # This a random walk (Wichman 2017 implementation) with path distance 
  # controlled by the two-parameter PCM friction model (Perla et al 1980)
  # for a sink-filled DEM.
  
  #__Terrain and Source Data__
  # dem:        The digital elevation model (SpatRast object)
  # xy:         The x,y coordinates of start (source) point in same CRS as
  #             the dem (Vector object; e.g. c(x,y))
  
  # __Random Walk Parameters__
  # slp_thresh: A slope threshold (∘) defining where divergent flow is allowed
  # exp_div:    The exponent of divergence that controls the amount of divergence 
  #             or lateral spreading in areas below the slope threshold
  # per_fct:    A persistence factor that controls the direction of movement 
  # walks:       Number of simulations/iterations of random walk from xy
  
  # __PCM Friction Model Parameters__
  # mu:         The sliding friction coefficient - single numeric value or 
  #             'SpatRaster' object with same resolution and extent as dem.
  # md:         The mass-to-drag ratio
  # int_vel:    Initial velocity (m/s)
  
  #---Random Walk Initialization---
  
  # Get cell (number) from the DEM raster of the source point location 
  cntr_cell <- terra::cellFromXY(dem, xy = xy)
  sim_paths <- vector(mode = "list", length = walks) # for storing each simulated RW path
  sim_velocity <- vector(mode = "list", length = walks) # for storing each simulated velocity
  # Get resolution of the DEM raster - for cell distance calculations
  r = terra::res(dem)
  
  # Check if mu is single value or spatial varying (raster)
  is_sp_mu <- methods::is(mu, "SpatRaster")
  
  # Get row and column of the start cell in the DEM raster
  rw <- terra::rowFromY(dem, xy[2])
  cl <- terra::colFromX(dem, xy[1])
  rowcol <- c(rw, cl)
  
  # Get the row columns of neighboring cells to the source point (xy)
  d <- adjCells(r, xy)
  ngh_cells <- terra::cellFromXY(dem, d)
  
  # Get the coordinates of the neighboring cells
  ngh_points <- terra::xyFromCell(dem, ngh_cells) # could make row/col instead to keep as integer
  cntr_point <- terra::xyFromCell(dem, cntr_cell)
  
  # Calculate the distance to neighboring cells (will always be the same for a given DEM raster)
  cell_dist <- apply(ngh_points, 1, euclideanDistance, p2 = cntr_point) # this is always going to be constant - should be done outside of loop for multiple
  
  # Convert the DEM raster to a matrix to help make the calculations quicker
  m_dem <- terra::as.matrix(dem, wide=TRUE) # !pull out of function or keep out of loop for multiple source points...
  
  if(is_sp_mu){
    mu <- terra::as.matrix(mu, wide=TRUE) 
  }
  
  #---Random Walk Iterations---
  # Loop through k repetitions/simulations (walks) creating a unique random walk path for each one
  
  for(k in 1:walks){
    
    #print(k)
    path_cells <- list()
    cntr_cell <- rowcol
    
    # for storing max velocity of PCM
    vel_cells <- list()
    
    #path_cells[[1]] <- cntr_cell # added to make sure start cell is in path cells
    path_cells[[1]] <- NA #
    
    i = 0 # i identifies step within walk
    n = 0
    prv_pos <- 9999
    v_p = int_vel # Initial velocity (m/s) for PCM model
    theta_p = 1
    
    
    while(n == 0){
      
      i = i + 1 # updates step within walk (iterative)
      
      d <- adjRowCol(cntr_cell) # updates cell neighbors for current step
      
      # below are breaks to stop the walk if step is at an edge of the raster
      if(cntr_cell[1]+1 > nrow(m_dem)){
        #vel_cells[[i]] <- 0
        break # stop if at edge of DEM raster
      }
      
      if(cntr_cell[2]+1 > ncol(m_dem)){
        #vel_cells[[i]] <- 0
        break # stop if at edge
      }
      
      # get elevation values for cell neighbors
      elv_values <- m_dem[d]
      
      if(any(is.na(elv_values))){
        #vel_cells[[i]] <- 0
        break # stop if at edge
      }
      
      if(length(elv_values) < 8){
        #vel_cells[[i]] <- 0
        break # stop if at edge
      }
      
      # get elevations of neighboring cells
      elv_ngh <- elv_values[1:8]
      cell_id <- 1:8
      elv_cntr <- m_dem[cntr_cell[1], cntr_cell[2]]
      
      # check if any of neighboring cells have lower elevation
      lower_elv <- elv_ngh < elv_cntr
      
      if(!any(lower_elv)){
        #vel_cells[[i]] <- 0
        break # stop if no lower elevations
      }
      
      # determine slope to neighbor cell
      beta_ngh <-  atan( (elv_cntr - elv_ngh)  / cell_dist) * 180/pi
      
      if(anyNA(beta_ngh)){
        #vel_cells[[i]] <- 0
        break # stop approaching NA value in DEM raster
        #^ can comment out if want to allow random walk with some NA values in
        # neighbor cells
      }
      
      # start by giving equal weighting factors (f) to each neighbor cell
      f <- rep(1, 8)
      
      # use persistence factor (per_fct) weighting factor if the local flow 
      # direction equals the previous flow direction
      if(prv_pos < 9){ 
        f[prv_pos] <- per_fct
      }
      
      # filter weighing factors by neighboring cells lower than local cell
      f <- f[lower_elv]
      cells <- cell_id[lower_elv]
      
      # get slope angle of neighboring cells lower than local cell
      beta_ngh <- beta_ngh[lower_elv]
      
      # determine slope value (gamma_i) based on the slope threshold (slp_thres)
      gamma_i <- tan(beta_ngh*pi/180) / tan(slp_thresh*pi/180)
      
      # determine probability (prob) for each neighbor cell to be randomly selected
      fj <- f * tan(beta_ngh*pi/180)
      prob <-  f*tan(beta_ngh*pi/180) / sum(fj, na.rm = TRUE)
      prob <- prob/sum(prob, na.rm = TRUE)
      
      # determine how close the slope to the steepest neighbor is to the 
      # slope threshold
      gamma_max <- max(gamma_i, na.rm = TRUE)
      
      if(gamma_max > 1){
        # if gamma > select only steepest neighbor - can have ties
        N <- cells[gamma_max == gamma_i]
        
        # break ties with random selection
        if(length(N) > 1){
          nxt_cell <- sample(N, size = 1)
          
        } else {
          nxt_cell <- N
        }
        
      } else {
        # otherwise multiple flow directions for debris flows (mfdf) criterion
        # (Gamma 2000)
        
        # neighbor cells are filtered by allowable flow spread by the exponent
        # of divergence (exp_div)
        N <- cells[gamma_i >= gamma_max^exp_div]
        
        # assign transition probabilities to this selection 
        trans_prob <- prob[gamma_i >= gamma_max^exp_div]
        
        # use weighted sampling
        if(length(N) > 1){
          nxt_cell <- sample(N, size = 1, prob = trans_prob)
          
        } else {
          nxt_cell <- N
        }
        
      }
      
      # store selected cell for reference as previous cell
      prv_pos <- nxt_cell
      
      # Apply PCM two-parameter friction model to determine when to stop
      # the random walk path
      
      # Change to spatial varying friction coefficient if mu is supplied as a raster
      if(is_sp_mu){
        mu_in <-  mu[d[nxt_cell,][1], d[nxt_cell,][2]]
      } else {
        mu_in <- mu
      }
      
      # calculate walk velocity for this instance
      v_i <- pcm(mu = mu_in, md = md, v_p = v_p, theta_p = theta_p, theta_i = beta_ngh[nxt_cell == cells], l = cell_dist[nxt_cell])
      
      # record velocity
      vel_cells[[i]] <- v_i
      
      # assign current velocity as previous velocity (for next step in the walk)
      v_p <- v_i
     
      # assign current slope angle as previous slope angle (for next step in the walk)
      theta_p = beta_ngh[nxt_cell == cells] 
      
      # store row column of next cell for the walk (to reference DEM location)
      path_cells[[i]] <- d[nxt_cell,]
      
      # define local cell for next step in the walk
      cntr_cell <- d[nxt_cell,]
      
      # stop if no velocity calculated
      
      if(is.nan(v_i)){
        # record velocity
        vel_cells[[i]] <- 0
        break
      }
      
    }
    
    # Store walk path for k-th iteration/simulation # removed b/c of memory issues...
    # path_row_col <- matrix(unlist(path_cells), ncol = 2, byrow = TRUE)
    # if(is.numeric(path_row_col)){
    #   sim_paths[[k]] <- apply(path_row_col, 1, function(x) terra::cellFromRowCol(dem, row = x[1], col = x[2]))
    # }
    
    sim_paths[[k]] <- matrix(unlist(path_cells), ncol = 2, byrow = TRUE)
    
    if(length(vel_cells)>0){
      # to deal with breaks with no velocity
      sim_velocity[[k]] <- round(unlist(vel_cells),3)
    }
    
    
  }
  
  # Remove walks that didn't go anywhere 
  
  sim_paths <- Filter(function(x) !all(is.na(x)), sim_paths)
  
  if(length(sim_paths)== 0){
    # if no paths found...
    walks_res <- list(
      start_cell =  terra::cellFromRowCol(dem, rowcol[1], rowcol[2]),
      cell_trav_freq = NA,
      cell_max_vel = NA,
      prob_connect = NA
    )
    
  } else {
    
    sim_velocity <- Filter(function(x) !all(is.na(x)), sim_velocity)
    
    all_indices <- unlist(lapply(sim_paths, function(path) {
      terra::cellFromRowCol(dem, path[, 1], path[, 2])
    }))
    
    vel_data <- data.frame(cell = all_indices, velocity = unlist(sim_velocity))
    
    # Aggregate to find the max velocity per cell
    max_velocity_by_cell <- aggregate(velocity ~ cell, data = vel_data, FUN = max)
    
    
    # Count occurrences of each cell index [freq of each cell traversed by cell size]
    cell_counts <- table(all_indices) 
    
    if(source_connect == TRUE){
      
      prob_connect <- sourceConnect(sim_paths = sim_paths, feature_mask = connect_feature, trials = walks)
      
      walks_res <- list(
        start_cell =  terra::cellFromRowCol(dem, rowcol[1], rowcol[2]),
        cell_trav_freq = cell_counts,
        cell_max_vel = max_velocity_by_cell$velocity,
        prob_connect = prob_connect
      )
      
      
    } else {
      
      walks_res <- list(
        start_cell =  terra::cellFromRowCol(dem, rowcol[1], rowcol[2]),
        cell_trav_freq = cell_counts,
        cell_max_vel = max_velocity_by_cell$velocity,
        prob_connect = NULL
      )
      
    }
    
  }
  

  
  return(walks_res)
  
}





makeSourceList <- function(source_xy)
{
  lapply(seq_len(nrow(source_xy)), function(i) matrix(source_xy[i, ], ncol = 2))
}

