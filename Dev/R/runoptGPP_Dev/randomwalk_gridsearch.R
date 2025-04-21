#' Grid search optimization for random walk runout path
#'
#' Computes the area under the receiver operating characteristic curve (AUROC) for
#'      runout paths simuluated using the random walk model component of the
#'      GPP tool in SAGA-GIS. The AUROC compares a runout polygon to the simulated
#'      path.
#' @param dem A DEM as a RasterLayer object
#' @param slide_plys Runout tracks as a SpatialPolygonsDataFrame
#' @param slide_src SSource points as a SpatialPointsDataFrame or source areas
#'      as a SpatialPolygonsDataFrame
#' @param slide_id Selects a single runout polygon from slide_plys by row
#' @param slp_v A vector of random walk slope thresholds - below lateral spreading is modeled
#' @param ex_v A vector or random walk exponents controlling lateral spread
#' @param per_v A vector or random walk persistence factors to weight flow direction consistency
#' @param gpp_iter Number of model iterations
#' @param buffer_ext (Optional) Defines buffer distance (in meters) around runout polygon
#'      to crop source DEM. This helps to reduce computational time
#' @param buffer_source (Optional) Can define a buffer distance (in meters) to extend source
#'      point to a source area
#' @param save_res (logical), if TRUE, will save results in current working directory
#' @param plot_eval If TRUE, will plot random walk path and runout polygon
#' @param saga_lib The initiated SAGA-GIS geoprocessor object
#' @return the area under the receiver operating characteristic
#' @details Runout is either simulated from a single source point or a buffered
#'      area round the source point.
#' @examples


rwGridsearch <- function(dem, slide_plys, slide_src,
                      slide_id = NULL, slp_v, ex_v, per_v,
                      gpp_iter = 1000, buffer_ext = 500, buffer_source = NULL,
                      save_res = FALSE, plot_eval = FALSE, saga_lib = NULL)

{
  
  # Coerce to spatial "sp" object
  if(class(slide_plys)[1] == "sf"){
    slide_plys = sf::as_Spatial(slide_plys)
  }
  
  if(class(slide_src)[1] == "sf"){
    slide_src = sf::as_Spatial(slide_src)
  }
  
  if(is.null(slide_id)){
    slide_id = 1
  }

  roc_result.nm <- paste("result_rw_roc_", slide_id, ".Rd", sep="")

  column.names <- ex_v
  row.names <- slp_v
  matrix.names <- per_v

  roc_result <- array(NA ,dim = c(length(slp_v), length(ex_v), length(per_v)),dimnames = list(row.names, column.names, matrix.names))

  #roc[row, col, matrix]
  #roc[rwslp, rwexp, rwper]

  for(j in 1:length(per_v)){
    for(i in 1:length(ex_v)){
      for(k in 1:length(slp_v)){
        #res[k, i, j] <- paste(pcmmd[k], rwexp[i], rwper[j] )
        roc <- rwPerformance(dem,
                             slide_plys = slide_plys,
                             slide_src = slide_src,
                             slide_id = slide_id,
                             slp = slp_v[k],
                             ex = ex_v[i],
                             per = per_v[j],
                             gpp_iter = gpp_iter,
                             buffer_ext = buffer_ext,
                             buffer_source = buffer_source,
                             plot_eval = plot_eval,
                             saga_lib = saga_lib)

        roc_result[k, i, j] <- roc

      }}}


  if(save_res){
    save(roc_result, file=roc_result.nm)
  }

  return(roc_result)
}


#' Get random walk grid search optimal parameters
#'
#' @param x A list of random walk grid search values
#' @param measure A measure (e.g. median, mean) to find optimal parameter across grid search space
#' @param from_save (Logical) if TRUE, will load save files from current working directory
#' @return A dataframe  with the optimal parameter set and AUROC performance

rwGetOpt <- function(x, measure = median, from_save = FALSE){

  if(from_save){

    x <- list()

    files <- list.files(pattern = "result_rw_roc")

    for(i in 1:length(files)){
      res_nm <- paste("result_rw_roc_", i, ".Rd", sep="")
      res_obj_nm <- load(res_nm)
      roc_result <- get(res_obj_nm)
      x[[i]] <- roc_result
    }


  }

  rwslp_vec <- as.numeric(dimnames(x[[1]])[[1]])
  rwexp_vec <- as.numeric(dimnames(x[[1]])[[2]])
  rwper_vec <- as.numeric(dimnames(x[[1]])[[3]])

  roc_list <- list()

  for(PER in 1:length(rwper_vec)){

    per_list <- list()

    for(i in 1:length(x)){
      res <- x[[i]]
      per_list[[i]] <- as.matrix(res[,,PER])

    }

    Y <- do.call(cbind, per_list)
    Y <- array(Y, dim=c(dim(per_list[[1]]), length(per_list)))

    Y.measure <- apply(Y, c(1, 2), measure, na.rm = TRUE)

    roc_list[[PER]] <- Y.measure

  }

  roc <- do.call(cbind, roc_list)
  roc <- array(roc, dim=c(dim(roc_list[[1]]), length(roc_list)))

  # Find which combination of RW parameters had the highest auroc
  ROCwh <- which(roc==max(roc), arr.ind=T)
  ROCwh
  roc[ROCwh]

  opt_params <- data.frame(
    rw_slp_opt = rwslp_vec[ROCwh[1]],
    rw_exp_opt = rwexp_vec[ROCwh[2]],
    rw_per_opt = rwper_vec[ROCwh[3]],
    rw_auroc = roc[ROCwh]
  )

  return(opt_params)

}




#' Get random walk runout distance grid search optimal parameters for single event
#'
#' @param x the results of the grid search
#' @return A dataframe with the optimal parameter set and performance

rwGetOpt_single <- function(x){

  slp_vec <- as.numeric(dimnames(x)[[1]])
  ex_vec <- as.numeric(dimnames(x)[[2]])
  per_vec <- as.numeric(dimnames(x)[[3]])

  wh_opt <- which(x==max(x), arr.ind = TRUE)

  if(nrow(wh_opt) > 1){
    wh_opt <- wh_opt[1,]
  }

  opt_param <- data.frame(
    rw_slp_opt = slp_vec[wh_opt][1],
    rw_exp_opt = ex_vec[wh_opt][2],
    rw_per_opt = per_vec[wh_opt][3],
    rw_auroc = x[wh_opt][1]
  )

  return(opt_param)

}

