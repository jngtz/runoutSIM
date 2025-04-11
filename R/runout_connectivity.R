sourceConnect <- function(sim_paths, feature_mask, trials = NULL) {
  # Handle empty sim_paths
  if (length(sim_paths) == 0) {
    return(0)  # No connectivity if no paths
  }
  
  # Determine the total number of paths to process
  N_total <- if (is.null(trials)) length(sim_paths) else min(trials, length(sim_paths))
  
  # Counter for the number of paths that intersect the feature
  N_intersect <- 0
  
  # Loop over each path
  for (i in 1:N_total) {
    path <- sim_paths[[i]]  # Get the current path
    
    # Check if the path intersects the feature (assuming feature_mask is a matrix)
    if (any(feature_mask[path[, 1], path[, 2]] == 1)) {
      N_intersect <- N_intersect + 1
    }
  }
  
  # Calculate the connectivity probability
  P_connect <- N_intersect / N_total
  
  return(P_connect)
}


# Create a feature for connectivity (object) analysis using the river
makeConnFeature <- function(x,y){
  object = terra::rasterize(vect(x), y)
  
  # Create feature mask outside of function (just raster converted to matrix)
  feature_mask <- as.matrix(object, wide=TRUE)
  feature_mask[is.na(feature_mask)] <- 0
  return(feature_mask)
}