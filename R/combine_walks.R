
walksToRaster <- function(x, dem, weights = NULL){
  #Takes output of runoutSim and makes into a grid of
  # traverse frequencies for each cell. Can be adjusted by weights
  # (vector of the same length as source cells used for runoutSim)
  
  runout_raster <- terra::rast(dem)
  terra::values(runout_raster) <- NA
  
  if(length(x) > 3){
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
