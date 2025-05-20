#' Convert Random Walk Results to Raster
#'
#' Aggregates output from `runoutSim()` into a raster showing how often each
#' cell was traversed during runout simulations. Optionally weights each walk
#' by a supplied vector (e.g., probability of initiation).
#'
#' @param x A list of runout simulation outputs (e.g., from `runoutSim()`).
#' @param dem A `SpatRaster` object representing the reference DEM.
#' @param method Method of combining overlaying runout paths: "freq" and "cdf_prob" calculate the traverse frequency 
#' and cumulative distribution function probabilities for simulations from a single source cell. "max_cdf_prob" and "avg_cdf_prob"
#' calculate the maximum/average empirical probabilites from CDF's applied to walks from invididual source cells.
#' @param weights Optional numeric vector of weights (same length as `x`) to scale traverse frequencies.
#'
#' @return A `SpatRaster` with cell values representing weighted or unweighted traverse frequencies.
#'
#' @examples
#' \dontrun{
#' # Convert simulated paths to a raster of transition frequencies
#' paths_raster <- walksToRaster(sim_paths, dem)
#' plot(paths_raster)
#' }
#' @export

walksToRaster <- function(x, dem, method = "freq", weights = NULL){
  
  # Method = freq, avg_cdf_prob, and max_cdf_prob, cdf_prob (for single source)

  runout_raster <- terra::rast(dem)
  terra::values(runout_raster) <- NA
  
  if(!(is.list(x) && !is.null(names(x)))){
    
    #remove any with NA cell traversed
    x <- x[!sapply(x, function(el) any(is.na(as.vector(el$cell_trav_freq))))]
    
    # for multiple walks from difference source points
    trav_freq <- sapply(x, function(x) x$cell_trav_freq)
    
    if(!is.null(weights)){
      
      # multiply weights (e.g. probability of being a source) by traverse frequency
      trav_freq <- lapply(seq_along(trav_freq), function(i) trav_freq[[i]] * weights[i])
      
    }
    
    
    if(method == "freq"){
      d <- unlist(trav_freq)

      # combine overlaying traverse frequencies
      combine_d <- tapply(d, names(d), sum, simplify = TRUE) 
  
    }
    
    if(method == "avg_cdf_prob"){
      
      trav_eprob <- sapply(trav_freq, function(x) cdfProb(x))
      
      d <- unlist(trav_eprob)
      
      # combine overlaying traverse frequencies
      combine_d <- tapply(d, names(unlist(trav_freq)), mean, simplify = TRUE) 
      
    }
    
    
    if(method == "max_cdf_prob"){
      trav_eprob <- sapply(trav_freq, function(x) cdfProb(x))
      
      d <- unlist(trav_eprob)
      
      # combine overlaying traverse frequencies
      combine_d <- tapply(d, names(unlist(trav_freq)), max, simplify = TRUE) 
      
    }
   
    
    runout_raster[as.numeric(names(combine_d))] <- as.numeric(combine_d)
    
  } else {
    # for walks from one source
    cell_counts <- x$cell_trav_freq
    
    if(!is.null(weights)){
      runout_raster <- runout_raster*weights
    }
    
    if(method == "cdf_prob"){
      cell_counts <- cdfProb(cell_counts)
    }
    
    # Assign counts to the raster
    runout_raster[as.numeric(names(x$cell_trav_freq))] <- as.numeric(cell_counts)
   
  }
  
  # Add name attibutes
  names(runout_raster) <- method
  varnames(runout_raster) <- "runout_walks"
  return(runout_raster)
  
}

#' Convert Connectivity Probabilities to Raster
#'
#' Creates a raster where each source cell is assigned the probability
#' of connecting to a downstream feature.
#'
#' @param x A list of outputs from multiple runout simulations, each containing
#' `start_cell` and `prob_connect`.
#' @param dem A `SpatRaster` (typically the DEM used in simulation) used to define the output raster grid.
#'
#' @return A `SpatRaster` with probability of connectivity assigned to each source cell.
#'
#' @examples
#' \dontrun{
#' # Convert connectivity results to a raster
#' conn_prob <- connToRaster(sim_paths, dem)
#' plot(conn_prob)
#' }
#' @export

connToRaster <- function(x, dem){
  
  #remove any with NA cell traversed
  x <- x[!sapply(x, function(el) any(is.na(as.vector(el$cell_trav_freq))))]
  
  prob_connect <- round(sapply(x, function(x) x$prob_connect),3)
  cell_index <- sapply(x, function(x) x$start_cell)
  
  
  conn_r <- terra::rast(dem)
  terra::values(conn_r) <- NA
  
  # Assign connectivity to cells
  conn_r[cell_index] <- prob_connect
  names(conn_r) <- "connectivity_prob"
  varnames(conn_r) <- "connectivity_prob"
  return(conn_r)
}

#' Convert Runout Velocities to Raster
#'
#' Aggregates maximum velocities from runout simulation paths into a raster,
#' using a specified summary method if multiple simulations overlap.
#'
#' @param x A list of simulation outputs containing `cell_max_vel` and `cell_trav_freq`.
#' @param dem A `SpatRaster` (DEM) used as reference.
#' @param method Summary function to apply across overlapping simulations (e.g., `"max"`, `"mean"`).
#'
#' @return A `SpatRaster` with cell values representing maximum velocity (in m/s).
#'
#' @examples
#' \dontrun{
#' # Convert velocity results to a raster
#' trav_vel <- velocityToRaster(sim_paths, dem)
#' plot(trav_vel)
#' }
#' @export

velocityToRaster <- function(x, dem, method = "max"){
  
  #remove any with NA cell traversed
  x <- x[!sapply(x, function(el) any(is.na(as.vector(el$cell_trav_freq))))]
  
  r <- terra::rast(dem)
  terra::values(r) <- NA
  
  if(!(is.list(x) && !is.null(names(x)))){
    # for multiple walks from difference source points
    cell_velocities <- sapply(x, function(x) x$cell_max_vel)
    cell_indicies <- sapply(x, function(x) as.numeric(names(x$cell_trav_freq)))
    
    vels <- unlist(cell_velocities)
    cells <- unlist(cell_indicies)
    
    # combine overlaying traverse frequencies
    combine_d <- tapply(vels, cells, get(method), simplify = TRUE) 
    
    
    r[as.numeric(names(combine_d))] <- as.numeric(combine_d)
    
  } else {

    cell_index <- as.numeric(names(x$cell_trav_freq))
    cell_velocity <- x$cell_max_vel
    
    r[cell_index] <- cell_velocity

  }
  
  # Add name attibutes
  names(r) <- "max_velocity_ms"
  varnames(r) <- method
  return(r)
  
  
}

#' Convert Raster Values to Cumulative Distribution
#'
#' Applies an empirical cumulative distribution function (ECDF) to raster values,
#' returning a raster where each cell reflects its percentile rank.
#'
#' @param x A `SpatRaster` object (e.g., from `walksToRaster` or `velocityToRaster`).
#'
#' @return A `SpatRaster` of the same extent and resolution with values between 0 and 1.
#'
#' @examples
#' \dontrun{
#' percentile_raster <- rasterCdf(runout_raster)
#' }
#' @export

rasterCdf <- function(x){
  
  x_ecdf <- stats::ecdf(terra::values(x))
  prob_x <- terra::setValues(x, x_ecdf(terra::values(x)))
  prob_x
  
}

#' Convert Raster Values to Cumulative Distribution
#'
#' Applies an empirical cumulative distribution function (ECDF) to raster values,
#' returning a raster where each cell reflects its percentile rank.
#'
#' @param x A vector
#'
#' @return A vector of empirical probabilities
#'
#' @examples
#' \dontrun{
#' test <- cdfProb(c(1,2,3,4,5))
#' }
#' @export

cdfProb <- function(x){
  x_ecdf <- stats::ecdf(x)
  prob_den <- x_ecdf(x)
  return(prob_den)
}

