
#' Random walk runout path performance
#'
#' Computes the area under the receiver operating characteristic curve (AUROC) for
#'      runout paths simuluated using the random walk model component of the
#'      GPP tool in SAGA-GIS. The AUROC compares a runout polygon to the simulated
#'      path.
#' @param dem A DEM as a RasterLayer object
#' @param slide_plys Runout tracks as a SpatialPolygonsDataFrame
#' @param slide_src Source points as a SpatialPointsDataFrame or source areas
#'      as a SpatialPolygonsDataFrame
#' @param slide_id Selects a single runout polygon from slide_plys by row to run GPP model
#' @param slp Random walk slope threshold - below lasteral spreading is modelled
#' @param ex Random walk exponent controlling lateral spread
#' @param per Random walk persistence factor to weight flow direction consistency
#' @param gpp_iter Number of model iterations
#' @param buffer_ext (Optional) Defines buffer distance (in meters) around runout polygon
#'      to crop source DEM. This helps to reduce computational time
#' @param buffer_source (Optional) Can define a buffer distance (in meters) to extend source
#'      point to a source area
#' @param plot_eval logical. If TRUE will plot random walk path and runout polygon
#' @param saga_lib The initiated SAGA-GIS geoprocessor object
#' @return The area under the receiver operating characteristic
#' @details Runout source can be either point or area.
#' @examples
#' \dontrun{
#' # Initialize a saga object
#' saga <- Rsagacmd::saga_gis()
#'
#' # Load elevation model (DEM)
#' dem <- raster::raster(system.file("extdata/elev_12_5m.tif", package="runout.opt"))
#'
#' # Load runout polygons and source points
#' runout_plys <- rgdal::readOGR(system.file("extdata/dflow_runout_ply.shp", package="runout.opt"))
#' source_pnts <- rgdal::readOGR(system.file("extdata/dflow_source_pnt.shp", package="runout.opt"))
#'
#' # Run GPP random walk model for a rounout polygon
#' rw <- rwPerformance(dem, slide_plys = runout_plys[1,], slide_src = source_pnts,
#'     slp = 30, ex = 3, per = 2,
#'     gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
#'     plot_eval = TRUE, saga_lib = saga)
#'
#' rw # returns AUROC
#'
#' }


rwPerformance <- function(dem, slide_plys, slide_src,
                          slide_id = 2, slp = 33, ex = 3, per = 2,
                          gpp_iter = 1000, buffer_ext = 500, buffer_source = NULL,
                          plot_eval = FALSE, saga_lib = NULL)
{

  # Coerce to spatial "sp" object
  if(class(slide_plys)[1] == "sf"){
    slide_plys = sf::as_Spatial(slide_plys)
  }
  
  if(class(slide_src)[1] == "sf"){
    slide_src = sf::as_Spatial(slide_src)
  }
  
  # If single runout polygon as input, assign slide_id value of 1
  if(length(slide_plys) == 1){
    slide_id <- 1
  }

  slide_plys$objectid <- 1:length(slide_plys)
  # Subset a single slide polygon
  slide_poly_single <- slide_plys[slide_id,]

  # Crop dem to slide polygon
  dem_grid <- raster::crop(dem, raster::extent(slide_poly_single) + buffer_ext)


  if(class(slide_src) == "SpatialPointsDataFrame"){
    # Subset corresponding source/start point of runout
    sel_over_start_point  <- sp::over(slide_src, slide_poly_single)
    sel_start_point <- slide_src[!is.na(sel_over_start_point$objectid),]

    if(!is.null(buffer_source)){
      # Create a buffer around source point to create a source/release area
      source_buffer <- sf::st_buffer(st_as_sf(sel_start_point), dist = buffer_source)
      source_grid <- raster::rasterize(source_buffer, dem_grid, field=1 )
      source_grid <- raster::mask(source_grid, slide_poly_single )
      source_plot <- raster::rasterToPolygons(source_grid)
    } else {
      # Just use source point
      source_plot <- sel_start_point
      source_grid <- raster::rasterize(matrix(sp::coordinates(sel_start_point)[1:2], ncol = 2), dem_grid, field = 1)
    }
  }

  if(class(slide_src) == "SpatialPolygonsDataFrame" ){
    sel_over_start_poly <- sp::over(slide_src, slide_poly_single)
    sel_start_poly <- slide_src[!is.na(sel_over_start_poly$objectid),]
    source_plot <- sel_start_poly
    source_grid <- raster::rasterize(sel_start_poly, dem_grid, field=1 )

  }


  # Run runout model using Rsagacmd (faster than RSAGA)
  if(is.null(saga_lib)){
    
    dem_terra <- terra::rast(dem_grid)
    
    if(!length(slide_plys) == 1){
      
      # Many source points
      source_l <- raster::xyFromCell(source_grid, which(values(source_grid) == 1))
      source_l <- lapply(1:nrow(source_l), function(i) matrix(source_l[i, ], nrow = 1, ncol = 2))
      
      rw_paths <- lapply(source_l, function(x) {
        runoutSim(dem = dem_terra, xy = x, mu = 0.00001, md = 500, 
                  slp_thresh = slp, exp_div = ex, per_fct = per, walks = gpp_iter,
                  source_connect = FALSE)})
      
    } else {
      
      # One source point
      source_plot <- sel_start_point
      source_l <- sp::coordinates(slide_src)
      source_grid <- raster::rasterize(matrix(sp::coordinates(sel_start_point)[1:2], ncol = 2), dem_grid, field = 1)
      
      rw_paths = runoutSim(dem = dem_terra, xy=source_l, mu = 0.08, md = 140, 
                           slp_thresh = slp, exp_div = ex, per_fct = per, walks = gpp_iter,
                           source_connect = FALSE)
      
    }

    rw_grid <- walksToRaster(x = rw_paths, dem = dem_terra)
    #rw_grid <- rw_grid/gpp_iter
    #plot(rw_grid)
    rw_grid <- raster::raster(rw_grid)
    
    pred_values <- terra::values(rw_grid)
    pred_values[is.na(pred_values)] <- 0

    
  } else {
    
    gpp <- saga_lib$sim_geomorphology$gravitational_process_path_model(dem = dem_grid, release_areas = source_grid,
                                                                       process_path_model = 1,
                                                                       rw_slope_thres = slp,
                                                                       rw_exponent = ex,
                                                                       rw_persistence = per,
                                                                       gpp_iterations = gpp_iter)
    

    rw_grid <- rasterCdf(gpp$process_area)
    
    pred_values <- raster::getValues(rw_grid)
    pred_values[is.na(pred_values)] <- 0
    
    
  }
  
  # AUROC
  
  slide_area <- raster::rasterize(slide_poly_single, rw_grid, field=1, background = 0)
  
  pred_area <- ROCR::prediction(predictions = pred_values, labels = raster::getValues(slide_area))
  
  auroc_area <- ROCR::performance(pred_area, "auc")
  roc <- auroc_area@y.values[[1]]


  # Plot evaluation results
  if(plot_eval){
    raster::plot(rw_grid,
         main = paste("id", slide_id,
                      "auroc", round(roc, digits=3), "\n",
                      "Ex", ex, "Pr", per, "Slp", slp),
         cex.main = 0.7, cex.axis = 0.7, cex=0.7)
    sp::plot(slide_poly_single, add=TRUE)
    sp::plot(source_plot, add=TRUE)
  }

  return(roc)

}


