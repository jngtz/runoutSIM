# testing max freq approach to cacluating combine prob - tested and
# is same results as method = prob - and method is prob is faster
# will just leave out max_freq method...


multi.paths_to_raster <- function(x,y,method="prob", walks = NULL, weights = NULL){
  # requires terra::
  # x:        a list of simulated runout paths
  # y:        the dem (terra::) used to determine x
  # method:   "prob" returns traverse probabilities, 
  #           "freq" returns traverse frequency,
  #           "max freq" returns max frequency,
  # weights:  a vector of weights to adjust freq or pobabilities
  #           of the same length as x
  # walks:    number of iterations or walks performed from each source
  
  # Checks if number of walks given if method is probability
  if(method == "prob" & is.null(walks)){
    stop('Must provide number of walks for method = "prob"') 
  }
  
  # Create empty raster to store results
  sim <- y
  sim[] <- 0
  
  # Define initial values
  sim_values <- terra::values(sim)
  prev_values <- sim_values[] 
  
  if(method == "freq" || method == "max freq"){
    prev_values[] <- 0
  } else {
    prev_values[] <- 1 # to combine no traverse probabilities
  }
  
  # Loop to combine the path 
  for (i in 1:length(rw_l)) {
    
    # Reset simulation values to zero
    if(method == "freq" || method == "max freq"){
      sim_values[] <- 0 
    } else {
      sim_values[] <- 1 
    }
    
    # Process current simulation
    sim_i <- rw_l[[i]]
    
    # Flatten all path indices into a single vector
    all_indices <- unlist(lapply(sim_i, function(path) {
      terra::cellFromRowCol(sim, path[, 1], path[, 2])
    }))
    
    # Count occurrences of each cell index
    traverse_count <- table(all_indices)
    
    # Adjust by source probability
    if(!is.null(weights)){
      traverse_count <- traverse_count * weights[i]
    }
    
    
    if(method == "freq"){
      sim_values <- prev_values + sim_values
      sim_values[as.numeric(names(traverse_count))] <- traverse_count
    } else if (method == "max freq"){
      sim_values <- pmax(prev_values,sim_values)
      sim_values[as.numeric(names(traverse_count))] <- traverse_count
    } else {  
      
      # Convert cell counts to no traverse probability
      traverse_count = 1 - (traverse_count/walks)
      
      sim_values <- prev_values * sim_values
      sim_values[as.numeric(names(traverse_count))] <- traverse_count
    }
    
    
    prev_values <- sim_values
    
  }
  
  # Apply values to empty raster
  if(method == "freq" || method == "max freq"){
    sim[] <- sim_values
  } else {
    sim[] <- 1-sim_values # to transform back to prob of traverse
  }
  
  # Remove values of 0
  sim[sim == 0] <- NA
  
  return(sim)
  
}

