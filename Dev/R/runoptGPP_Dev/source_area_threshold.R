#' Raster binary reclassification
#'
#' Reclassifies a RasterLayer object based on range of values.
#' @param x A RasterLayer
#' @param value_range A vector with the min and max value for reclassification
#' @return A RasterLayer with two classes.

rasterThreshold <- function(x, value_range = c(0.9, Inf)){
  if(value_range[2] != Inf){

    m <- c(-Inf, value_range[1], NA, value_range, 1, value_range[2], Inf, NA)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)

  } else {

    m <- c(-Inf, value_range[1], NA, value_range, 1)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)

  }

  rc <- raster::reclassify(x, rclmat)
  rc

}


#' Runout process areas from source area prediction
#'
#' Thresholds a prediction raster of source areas (probabilities) and computes
#' corresponding runout process areas.
#' @param source_pred RasterLayer with probability of being runout source area
#' @param dem A DEM as a RasterLayer object
#' @param source_threshold A cutoff of value of source area prediction
#' @param rw_slp Random walk slope threshold - below lateral spreading is modeled
#' @param rw_exp Random walk exponent controlling lateral spread
#' @param rw_per Random walk persistence factor to weight flow direction consistency
#' @param pcm_mu PCM model sliding friction coefficient
#' @param pcm_md PCM model mass-to-drag ratio (m)
#' @param gpp_iter Model iterations
#' @param saga_lib The initiated SAGA-GIS geoprocessor object
#' @param file_name (optional) Can give a file name to save process area raster as GeoTiff
#' @return A RasterLayer with runout process areas

runoutPareaPredict <- function(source_pred, dem, source_threshold,
                               rw_slp, rw_exp, rw_per, pcm_mu, pcm_md,
                               gpp_iter = 1000, saga_lib, file_name = NULL){

  dem_grid <- dem

  source_class <- rasterThreshold(source_pred, value_range = c(source_threshold, Inf))

  #just in case no source cells are present...
  if(raster::freq(source_class, value = 1) != 0) {
    gpp <- saga_lib$sim_geomorphology$gravitational_process_path_model(dem = dem_grid,
                                                                   release_areas = source_class,
                                                                   process_path_model = 1,
                                                                   rw_slope_thres = rw_slp,
                                                                   rw_exponent = rw_exp,
                                                                   rw_persistence = rw_per,
                                                                   gpp_iterations = gpp_iter,
                                                                   friction_model = 5,
                                                                   friction_mu = pcm_mu,
                                                                   friction_mass_to_drag = pcm_md)

    if(!is.null(file_name)){
      write_name <- paste0(file_name, ".tif")
      raster::writeRaster(gpp$process_area, paste0(file_name), format = "GTiff", overwrite = TRUE)
    }
  }
  return(gpp$process_area)
}


#' Compute area under ROC curve for process areas
#'
#' @param process_area RasterLayer of GPP process area
#' @param slide_area Runout tracks as a RasterLayer
#' @param smp_mask A mask of sample areas as a RasterLayer
#' @param smp_size Size of random sample of runount and non-runout cells
#' @param seed Specific specific random seed
#' @return AUROC

rocParea <- function(process_area, slide_area, smp_mask, smp_size = 1000, seed = 1234){

  set.seed(seed)
  rescale_process_area <- rasterCdf(process_area)
  # AUROC


  #noslide_area <- raster::rasterize(slide_polys, dem, field=NA, background = 1)
  #noslide_area <- raster::mask(noslide_area, dem)
  noslide_area <- smp_mask
  names(noslide_area) <- "mask"

  set.seed(seed)
  noslide_smp <- raster::sampleRandom(noslide_area, smp_size, sp=TRUE)
  noslide_smp$act = 0
  noslide_smp$mask = NULL

  slide_smp <- raster::sampleRandom(slide_area, smp_size, sp=TRUE)
  slide_smp$layer = NULL
  slide_smp$act = 1

  smp <- rbind(slide_smp, noslide_smp)
  smp$pred <- raster::extract(rescale_process_area, smp)
  smp$pred[is.na(smp$pred)] <- 0

  pred_smp <- ROCR::prediction(predictions = smp$pred, labels = smp$act)
  #perf_smp <- performance(pred_smp, "tpr", "fpr")
  auroc_smp <- ROCR::performance(pred_smp, "auc")
  roc <- auroc_smp@y.values[[1]]

  return(roc)

}


#' Compute percent area covered by process areas
#'
#' @param process_area RasterLayer of GPP process area
#' @param dem A DEM as a RasterLayer
#' @return Area covered by process area

areaPerParea <- function(process_area, dem){
  v_dem <- raster::getValues(dem)
  n_dem_cells <- length( v_dem[!is.na(v_dem)] )

  v_all <- raster::getValues(process_area) # process area
  v_parea <- v_all[!is.na(v_all)]
  n_parea_cells <- length(v_parea)
  per_area <- n_parea_cells / n_dem_cells *100
  per_area
}

