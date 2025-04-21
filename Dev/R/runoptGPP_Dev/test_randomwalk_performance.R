

runout_ply = st_as_sf(test_polygons[1,])
runout_src = st_filter(st_as_sf(source_points), st_as_sf(runout_ply))

dem = dem
runout_id = 1
slp = 33
ex = 3
per = 2
gpp_iter = 1000
buffer_ext = 500
buffer_source = NULL
plot_eval = FALSE

rwPerformance <- function(dem, runout_ply, runout_src,
                          runout_id = 1, slp = 33, ex = 3, per = 2,
                          walks = 1000, buffer_ext = 500, buffer_source = NULL,
                          plot_eval = FALSE)
{
  
  
  # Crop dem to runout polygon for faster computation
  dem_crop <- terra::crop(dem, terra::ext(runout_ply) + buffer_ext)
  start_point <- sf::st_coordinates(runout_src)
  
  # Run runout model using Rsagacmd (faster than RSAGA)
  
  rw_paths <- pcmRW(dem = dem_crop, start_point, int_vel = 1, 
                    slp_thresh = slp, exp_div = ex, per_fct = per,
                    mu = 0.0001, md = 0.001, reps = walks)
  
  rw_grid <- paths_to_raster(dem_crop, rw_paths)
  rw_grid <- rw_grid/walks
  plot(rw_grid)
  
  pred_values <- terra::values(rw_grid)
  pred_values[is.na(pred_values)] <- 0
  
  obs_area <- terra::rasterize(runout_ply, rw_grid, field=1, background = 0)
  obs_values <- terra::values(obs_area)
  
  pred_area <- ROCR::prediction(predictions = pred_values, labels = obs_values)
  
  auroc_area <- ROCR::performance(pred_area, "auc")
  roc <- auroc_area@y.values[[1]]
  
  # Plot evaluation results
  if(plot_eval){
    terra::plot(rw_grid,
                 main = paste("id", runout_id, "auroc", round(roc, digits=3), "\n",
                              "Ex", ex, "Pr", per, "Slp", slp),
                 cex.main = 0.7, cex.axis = 0.7, cex=0.7)
    plot(st_geometry(runout_ply), add=TRUE)
    plot(st_geometry(runout_src), add=TRUE)
  }
  
  return(roc)
  
}


rw.roc <- rwPerformance(dem = dem, runout_ply = runout_ply, runout_src = runout_src,
              runout_id = "test", slp = 40, ex = 3, per = 1,
              walks = 1000, buffer_ext = 500, buffer_source = NULL,
              plot_eval = TRUE)







rwGridsearch <- function(dem, runout_ply, runout_src,
                         slide_id = NULL, slp_v, ex_v, per_v,
                         gpp_iter = 1000, buffer_ext = 500, buffer_source = NULL,
                         save_res = FALSE, plot_eval = FALSE)
  
{
  
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
        roc <- rwPerformance(dem = dem, 
                             runout_ply = runout_ply, 
                             runout_src = runout_src,
                             runout_id = "test", 
                             slp = slp_v[k],
                             ex = ex_v[i],
                             per = per_v[j],
                             walks = gpp_iter, 
                             buffer_ext = buffer_ext, 
                             buffer_source = buffer_source,
                             plot_eval = plot_eval)
          
        roc_result[k, i, j] <- roc
        
      }}}
  
  
  if(save_res){
    save(roc_result, file=roc_result.nm)
  }
  
  return(roc_result)
}

steps <- 3
rwexp_vec <- seq(1.3, 3, len=steps)
rwper_vec <- seq(1.5, 2, len=steps)
rwslp_vec <- seq(20, 40, len=steps)

rw_gridsearch <- rwGridsearch(dem, runout_ply = runout_ply, runout_src = runout_src,
             slide_id = 1, slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
             gpp_iter = 1000, buffer_ext = 500, buffer_source = 50, save_res = TRUE,
             plot_eval = TRUE)


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

rwGetOpt_single(rw_gridsearch)




# PCM ##########################################################################

dem
runout_ply
runout_src 
slide_id = 1
rw_slp = 33
rw_ex = 3
rw_per = 2
pcm_mu = 0.3
pcm_md = 75
buffer_ext = 500
buffer_source = 50
gpp_iter = 1000
predict_threshold = 0.5
plot_eval = FALSE
return_features = FALSE


errMinBboxLength <- function(obs_poly, pred_raster, dem){
  #Calculates the relative error b/w bbox length estimates of slides
  sp.pred <- raster::rasterToPolygons(pred_raster, n = 4, dissolve = TRUE, na.rm=TRUE)
  geom_pred <- runoutGeom(sp.pred, elev = dem)
  geom_act <- runoutGeom(obs_poly, elev = dem)
  
  err <- geom_pred$length - geom_act$length
  rel_err <- abs(geom_act$length - geom_pred$length) / geom_act$length
  rel_diff<- (geom_pred$length - geom_act$length) / geom_act$length
  
  return(
    list(rel_error = rel_err,
         rel_difference = rel_diff,
         error = err
    ))
}


pcmPerformance <- function(dem, runout_ply, runout_src, slide_id = 1,
                           rw_slp = 33, rw_ex = 3, rw_per = 2,
                           pcm_mu = 0.3, pcm_md = 75,
                           buffer_ext = 500, buffer_source = 50, gpp_iter = 1000,
                           predict_threshold = 0.5, plot_eval = FALSE,
                           return_features = FALSE)
{
  
  # If single runout polygon as input, assign slide_id value of 1
  if(length(runout_ply) == 1){
    slide_id <- 1
  }
  
  runout_ply$objectid <- 1:length(runout_ply)
  # Subset a single slide polygon
  runout_src_coord <- st_coordinates(runout_src)
  
  # Crop dem to slide polygon
  dem_crop <- terra::crop(dem, terra::ext(runout_ply) + buffer_ext)
  
  
  # Run runout model using Rsagacmd (faster than RSAGA)
  gpp <- pcmRW(dem = dem_crop, runout_src_coord, int_vel = 1, 
                           slp_thresh = rw_slp, exp_div = rw_ex, per_fct = rw_per,
                           mu = pcm_mu, md = pcm_md, reps = gpp_iter)

  # remove na's
  gpp <- gpp[!sapply(gpp, function(x) any(is.na(x)))]
  
  rw_grid <- paths_to_raster(dem_crop, gpp)
  rw_grid <- rw_grid/gpp_iter
  plot(rw_grid)
  
  pred_values <- terra::values(rw_grid)
  pred_values[is.na(pred_values)] <- 0
  
  
  
  # Rescale to values from 0 to 1
  #rescale_process_area <- rasterCdf(gpp$process_area)
  
  # AUROC
  obs_area <- terra::rasterize(runout_ply, rw_grid, field=1, background = 0)
  obs_values <- terra::values(obs_area)
  
  pred_area <- ROCR::prediction(predictions = pred_values, labels = obs_values)
  
  auroc_area <- ROCR::performance(pred_area, "auc")
  roc <- auroc_area@y.values[[1]]
  # Length loss
  errMinBox <- errMinBboxLength(obs_poly = as_Spatial(runout_ply),
                                pred_raster = raster::raster(rw_grid),
                                dem = dem)
  
  length_relerr <- errMinBox[["rel_error"]]
  length_reldiff <- errMinBox[["rel_difference"]]
  length_error <- errMinBox[["error"]]
  
  # Plot evaluation results
  if(plot_eval){
    sp::plot(rescale_process_area,
             main = paste("id", slide_id,
                          "roc", round(roc, digits=2), "\n",
                          "err.L", round(length_relerr, digits = 2),
                          "Ex", rw_ex, "Pr", rw_per, "Mu", pcm_mu, "Md", pcm_md),
             cex.main = 0.7, cex.axis = 0.7, cex=0.7)
    sp::plot(slide_poly_single, add=TRUE)
    sp::plot(source_plot, add = TRUE)
  }
  
  
  if(return_features){
    return(
      list(id = slide_id,
           roc = roc,
           length.relerr = length_relerr,
           length.reldiff = length_reldiff,
           length.error = length_error,
           dem = dem_grid,
           actual.poly = slide_poly_single,
           source.pnt = sel_start_point,
           source.area = source_grid,
           gpp.parea = gpp$process_area,
           gpp.stop = gpp$stop_positions,
           gpp.maxvel = gpp$max_velocity))
  } else {
    return(
      list(id =  slide_id,
           roc = roc,
           length.relerr = length_relerr,
           length.reldiff = length_reldiff,
           length.error = length_error))
  }
  
}



