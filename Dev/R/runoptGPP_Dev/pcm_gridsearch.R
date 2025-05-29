#' Grid search optimization for PCM runout distance
#'
#' Computes performance measures for
#'      runout paths simuluated using the random walk and PCM model components of the
#'      GPP tool in SAGA-GIS.
#' @param dem A DEM as a RasterLayer object
#' @param slide_plys Runout tracks as a SpatialPolygonsDataFrame
#' @param slide_src Source points as a SpatialPointsDataFrame or source areas
#'      as a SpatialPolygonsDataFrame
#' @param slide_id Selects a single runout polygon from slide_plys by row
#' @param rw_slp Random walk slope threshold - below lateral spreading is modelled
#' @param rw_ex Random walk exponent controlling lateral spread
#' @param rw_per Random walk persistence factor to weight flow direction consistency
#' @param pcm_mu_v A vector of PCM model sliding friction coefficients
#' @param pcm_md_v A vector of PCM model mass-to-drag ratios (m)
#' @param gpp_iter Number of model iterations
#' @param buffer_ext (Optional) Defines buffer distance (in meters) around runout polygon
#'      to crop source DEM. This helps to reduce computational time when working
#'      with large regions.
#' @param buffer_source (Optional) Can define a buffer distance (in meters) to extend source
#'      point to a source area
#' @param predict_threshold A cutoff value to define what quantile of simulated runout
#'      frequencies is the predicted runout.
#' @param save_res (Logical) if TRUE, will save results in current working directory
#' @param plot_eval (Logical) ff TRUE, will plot simulated runout and runout polygon
#' @param saga_lib The initiated SAGA-GIS geoprocessor object
#' @return A list with performances for each parameter combinations across grid search space.
#'      The performances measures include relative error, relative difference, error and AUROC.

pcmGridsearch <- function(dem,
                       slide_plys, slide_src, slide_id = NULL,
                       rw_slp, rw_ex, rw_per,
                       pcm_mu_v, pcm_md_v,
                       gpp_iter = 1000,
                       buffer_ext = 500, buffer_source = NULL,
                       predict_threshold = 0.5,
                       save_res = FALSE, plot_eval = FALSE, saga_lib = NULL){

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

  column.names <- pcm_md_v
  row.names <- pcm_mu_v

  roc_result <- array(NA ,dim = c(length(pcm_mu_v),length(pcm_md_v)),dimnames = list(row.names, column.names))
  relerr_length_result <- roc_result
  reldiff_length_result <- roc_result
  err_length_result <- roc_result

  #roc[row, col]
  #roc[pcmmu, pcmmd]

  for(i in 1:length(pcm_md_v)){
    for(k in 1:length(pcm_mu_v)){

      result <- pcmPerformance(dem = dem,
                                 slide_plys = slide_plys,
                                 slide_src = slide_src,
                                 slide_id = slide_id,
                                 rw_slp = rw_slp,
                                 rw_ex = rw_ex,
                                 rw_per = rw_per,
                                 pcm_mu = pcm_mu_v[k],
                                 pcm_md = pcm_md_v[i],
                                 buffer_ext = buffer_ext,
                                 buffer_source = buffer_source,
                                 gpp_iter = gpp_iter,
                                 predict_threshold = predict_threshold,
                                 plot_eval = plot_eval,
                                 saga_lib = saga_lib)

      roc_result[k, i] <- result$roc
      relerr_length_result[k, i] <- result$length.relerr
      reldiff_length_result[k, i] <- result$length.reldiff
      err_length_result[k, i] <- result$length.error

      #roc_result[k, i] <- paste(pcmmd[i], pcmmu[k])
    }}

  result_pcm <- list(
    auroc = roc_result,
    relerr = relerr_length_result,
    reldiff = reldiff_length_result,
    error = err_length_result)

  if(save_res){
    save(result_pcm, file = paste("result_pcm_gridsearch_", slide_id, ".Rd", sep=""))
  }

  return(result_pcm)

}




#' Get PCM runout distance grid search optimal parameters
#'
#' @param x A list of PCM grid search values
#' @param performance Performance measure "relerr" relative error
#' @param measure A measure (e.g. "median", "mean") to find optimal parameter across grid search space
#' @param from_save (Logical) if TRUE, will load save files form current working directory
#' @param plot_opt (Logical) if TRUE, will plot performance across grid search space
#' @return A dataframe  with the optimal parameter set and performance

pcmGetOpt <- function(x, performance = "relerr", measure = "median",
                      from_save = FALSE, plot_opt = FALSE){

  if(from_save){
    x <- list()

    files <- list.files(pattern = "result_pcm_gridsearch_")

    for(i in 1:length(files)){
      res_nm <- paste("result_pcm_gridsearch_", i, ".Rd", sep="")
      res_obj_nm <- load(res_nm)
      result_pcm <- get(res_obj_nm)
      x[[i]] <- result_pcm
    }

  }

  pcm_md_vec <- as.numeric(colnames(x[[1]][[1]]))
  pcm_mu_vec <- as.numeric(rownames(x[[1]][[1]]))

  n_train <- length(x)

  err_list <- list()
  roc_list <- list()

  for(i in 1:n_train){

    err_list[[i]] <- x[[i]][[performance]]
    roc_list[[i]] <- x[[i]][['auroc']]

  }

  err <- apply(simplify2array(err_list), 1:2, get(measure))
  roc <- apply(simplify2array(roc_list), 1:2, get(measure))

  err_wh <- which(err==min(err), arr.ind=TRUE)
  err[err_wh]

  # Use AUROC for tie breaking
  if(length(err_wh) > 2){
    err_wh <- err_wh[which(roc[err_wh]==max(roc[err_wh]), arr.ind=TRUE),]
  }

  # If still no tie break, take first one...
  if(length(err_wh) > 2){
    err_wh <- err_wh[1,]
  }

  opt_md <- pcm_md_vec[err_wh[2]] #col
  opt_mu <- pcm_mu_vec[err_wh[1]] #row

  opt_gpp_par <- data.frame(
    pcm_mu = opt_mu,
    pcm_md = opt_md
    )

  opt_gpp_par[paste0(measure, "_", performance)] <- err[err_wh[1], err_wh[2]]
  opt_gpp_par[paste0(measure, "_", "auroc")] <- roc[err_wh[1], err_wh[2]]

  if(plot_opt){

    err_df <- reshape2::melt(err)

    gg <- ggplot2::ggplot(data = err_df, ggplot2::aes(x=err_df$Var2, y=err_df$Var1, z=err_df$value)) +
      ggplot2::geom_tile(ggplot2::aes(fill = err_df$value)) +

      ggplot2::ylab(expression(paste("Sliding friction coefficient"))) +
      ggplot2::xlab("Mass-to-drag ratio (m)") +

      ggplot2::labs(fill="Median relative\nrunout distance\nerror") +

      ggplot2::scale_fill_viridis_c(direction = 1) +
      ggplot2::theme_light() +
      ggplot2::theme(text = ggplot2::element_text(family = "Arial", size = 8),
                     axis.title = ggplot2::element_text(size = 9),
                     axis.text = ggplot2::element_text(size = 8))

   print(gg)
  }

  return(opt_gpp_par)

}


#' Get PCM runout distance grid search optimal parameters for single event
#'
#' @param x A list of PCM grid search values
#' @param performance Performance measure "relerr" relative error
#' @return A dataframe  with the optimal parameter set and performance

pcmGetOpt_single <- function(x, performance = "relerr"){

  pcm_md_vec <- as.numeric(colnames(x[[1]]))
  pcm_mu_vec <- as.numeric(rownames(x[[1]]))

  err <- x[[performance]]
  roc <- x[['auroc']]

  err_wh <- which(err==min(err), arr.ind=TRUE)
  err[err_wh]

  # Use AUROC for tie breaking
  if(length(err_wh) > 2){
    err_wh <- err_wh[which(roc[err_wh]==max(roc[err_wh]), arr.ind=TRUE),]
  }


  if(length(err_wh) > 2){
    err_wh <- err_wh[1,]
  }

  opt_md <- pcm_md_vec[err_wh[2]] #col
  opt_mu <- pcm_mu_vec[err_wh[1]] #row


  opt_gpp_par <- data.frame(
    pcm_mu = opt_mu,
    pcm_md = opt_md
  )

  opt_gpp_par[paste0(performance)] <- err[err_wh[1], err_wh[2]]
  opt_gpp_par$auroc <- roc[err_wh[1], err_wh[2]]

  opt_gpp_par
}


#' Get PCM runout distance grid search values
#'
#' @param x A list of PCM grid search values
#' @param performance Performance measure "relerr" relative error
#' @param measure A measure (e.g. "median", "mean") to aggregate model performances from all slides
#' @param from_save (Logical) if TRUE, will load save files form current working directory
#' @return A matrix of grid search space with aggregated performance values

pcmGetGrid <- function(x, performance = "relerr", measure = "median",
                      from_save = FALSE){

  if(from_save){
    x <- list()

    files <- list.files(pattern = "result_pcm_gridsearch_")

    for(i in 1:length(files)){
      res_nm <- paste("result_pcm_gridsearch_", i, ".Rd", sep="")
      res_obj_nm <- load(res_nm)
      result_pcm <- get(res_obj_nm)
      x[[i]] <- result_pcm
    }

  }

  pcm_md_vec <- as.numeric(colnames(x[[1]][[1]]))
  pcm_mu_vec <- as.numeric(rownames(x[[1]][[1]]))

  n_train <- length(x)

  err_list <- list()
  roc_list <- list()

  for(i in 1:n_train){

    err_list[[i]] <- x[[i]][[performance]]

  }

  err <- apply(simplify2array(err_list), 1:2, get(measure))

  err
}



