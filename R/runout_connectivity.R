#' Create Connectivity Feature Mask
#'
#' Converts a spatial feature (e.g., a polygon `sf` object) into a binary matrix
#' mask aligned with the reference DEM used for runout simulations.
#'
#' @param x An `sf` object representing the feature of interest (e.g., a polygon).
#' @param y A `SpatRaster` (e.g., DEM) used as the reference for rasterization.
#'
#' @return A binary matrix the same size as `y`, with `1` indicating presence of the feature and `0` elsewhere.
#'
#' @examples
#' feature_mask <- makeConnFeature(my_sf_polygon, my_dem)
makeConnFeature <- function(x,y){

  object = terra::rasterize(vect(x), y)
  
  # Create feature mask outside of function (just raster converted to matrix)
  feature_mask <- as.matrix(object, wide=TRUE)
  feature_mask[is.na(feature_mask)] <- 0
  return(feature_mask)
}


#' Estimate Probability of Source-Feature Connectivity
#'
#' Calculates the proportion of simulated paths (e.g., from random walk runout simulations)
#' that intersect a binary connectivity feature mask.
#'
#' @param sim_paths A list of matrices, each containing (row, col) indices of cells visited by a simulation path.
#' @param feature_mask A binary matrix indicating presence (`1`) of a connectivity feature (from `makeConnFeature()`).
#' @param trials Optional integer specifying how many simulation paths to consider. Defaults to all paths.
#'
#' @return A numeric value between 0 and 1 indicating the estimated probability of connection.
#'
#' @examples
#' p_connect <- sourceConnect(sim_paths, feature_mask, trials = 1000)
sourceConnect <- function(sim_paths, feature_mask, trials = NULL) {
  # This function works quickly by exploring for matches in (row cols)
  # between the connectivity feature and the walks. It counts the number
  # of walks in N total walks intersects with the feature mask (connected feature)
  
  # sim_paths, is list of row,col of cells that were traversed during the random walk
  # feature_mask is the prepared layer from makeConnFeature (a matrix indicating which
  # cells in the reference DEM the connectivity feature overlaps)
  
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


