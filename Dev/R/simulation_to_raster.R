

walksToRaster <- function(x, dem, weights = NULL){

  runout_raster <- terra::rast(dem)
  terra::values(runout_raster) <- NA
  
  if(length(x) > 4){
    # for multiple walks from difference source points
    trav_freq <- sapply(x, function(x) x$cell_trav_freq)
    
    if(!is.null(weights)){
      
      # multiply weights (e.g. probability of being a source) by traverse frequency
      trav_freq <- lapply(seq_along(trav_freq), function(i) trav_freq[[i]] * weights[i])
      
    }
    
    d <- unlist(trav_freq)
    
    # combine overlaying traverse frequencies
    combine_d <- tapply(d, names(d), sum, simplify = TRUE) 
    
    
    runout_raster[as.numeric(names(combine_d))] <- as.numeric(combine_d)
    
  } else {
    # for walks from one source
    cell_counts <- x$cell_trav_freq
    
    # Assign counts to the raster
    runout_raster[as.numeric(names(cell_counts))] <- as.numeric(cell_counts)
    
    if(!is.null(weights)){
      runout_raster <- runout_raster*weights
    }
    
  }
  
  # Add name attibutes
  names(runout_raster) <- "traverse_freq"
  varnames(runout_raster) <- "runout_walks"
  return(runout_raster)
  
}

connToRaster <- function(x, dem){
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


velocityToRaster <- function(x, dem, method = "max"){
  
  r <- terra::rast(dem)
  terra::values(r) <- NA
  
  if(length(x) > 4){
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


rasterCdf <- function(x){
  
  x_ecdf <- stats::ecdf(terra::values(x))
  prob_x <- terra::setValues(x, x_ecdf(terra::values(x)))
  prob_x
  
}
