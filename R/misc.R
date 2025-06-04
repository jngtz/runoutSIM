#' Extract DEM elevations for a given quantile within a polygon
#'
#' Identifies and masks the top or bottom quantile of elevation values 
#' within a given polygon overlaying a digital elevation model (DEM).
#'
#' @param dem A `SpatRaster` object representing a digital elevation model (DEM).
#' @param runout_ply A `SpatVector` polygon (or compatible `SpatialPolygons*` object) defining the area from which to extract elevation values.
#' @param quant A numeric value between 0 and 1 indicating the proportion (quantile) of elevation values to extract. Defaults to 0.05 (top/bottom 5%).
#' @param upper Logical. If `TRUE`, the function selects the highest elevation percentile; if `FALSE`, it selects the lowest. Defaults to `TRUE`.
#'
#' @return A `SpatRaster` object where selected cells are marked with 1 and all others are set to `NA`.
#'
#' @details This function is useful for identifying areas of extreme elevation within a polygon, such as potential landslide release zones (upper) or deposition zones (lower), based on percentile thresholds.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dem <- rast("dem.tif")
#' runout_ply <- vect("runout_area.shp")
#' release_zone <- elevPercentile(dem, runout_ply, per = 0.05, upper = TRUE)
#' }
#'
#' @export

elevQuantile <- function(dem , runout_ply, quant = 0.05, upper = TRUE){
  
  #per: percentile
  #dem: raster object elevation model
  #runout_ply: spatialPolygon
  
  slide_elev <- terra::extract(dem, runout_ply, cells = TRUE, df = TRUE)
  names(slide_elev) <- c("ID", "elev", "cell")
  
  num_cells <- nrow(slide_elev)
  cells_per <- quant * num_cells
  
  if(cells_per < 1){
    cells_per = 1
  }
  
  if(upper == TRUE){
    elev_ids <- as.numeric(rownames(slide_elev[order(slide_elev$elev, decreasing = TRUE),])[1:cells_per])
  } else {
    elev_ids <- as.numeric(rownames(slide_elev[order(slide_elev$elev, decreasing = FALSE),])[1:cells_per])
  }
  
  elev_cells <- slide_elev$cell[elev_ids]
  r_quantile <- terra::setValues(dem, NA)
  
  r_quantile[elev_cells] <- 1
  
  # Add name attibutes
  if(upper){
    sel_quant <- "upper"
  } else {
    sel_quant <- "lower"
  }
  
  names(r_quantile) <- paste("elev", quant, sel_quant, "quantile", sep = "_")
  varnames(r_quantile) <- paste("elev", quant, sel_quant, "quantile", sep = "_")
  return(r_quantile)
  
}