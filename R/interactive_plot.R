require(leaflet)
require(leafem)
require(viridis)
require(raster)
require(sf)
require(terra)

#' Plot Spatial Data Using Leaflet
#'
#' Creates an interactive `leaflet` map from raster or `sf` vector data, with options to customize basemaps, styling, and popups.
#'
#' This function is designed to handle both standalone and piped usage. If the first argument is a raster or `sf` object and `data` is `NULL`, it automatically reassigns the input appropriately. It can visualize raster values with color palettes and legends, or render vector data (points, lines, polygons) with informative popups.
#'
#' @param m Optional existing `leaflet` map object. If `NULL`, a new map is initialized.
#' @param data A `Raster*`, `terra::SpatRaster`, or `sf` object to be plotted. If not provided, `m` will be interpreted as the data.
#' @param group_layers Character vector of existing overlay groups. Used to maintain group-layer visibility toggles.
#' @param label A character string for the layer label and legend title. If `NULL`, will be auto-generated from the object name.
#' @param opacity Numeric (0–1) for layer transparency. Defaults to `0.5`.
#' @param color Color used for vector geometries (ignored for rasters). Defaults to `"black"`.
#' @param fill_color Color used for fill of vector geometries (ignored for rasters). Defaults to `color` parameter.
#' @param radius Numeric size of circle markers for point geometries. Defaults to `3`.
#' @param weight Line or border thickness for vector geometries. Defaults to `2`.
#' @param palette Color palette name used with `leaflet::colorNumeric()` for raster coloring. Defaults to `"viridis"`. If categorical values, supply a list - e.g.  list(classes = 1, colors = "#99d2ff", labels = "Water"))
#' @param basemaps Character vector of tile provider names (from `leaflet::providers`) to include as base layers. Defaults to `c("Esri.WorldImagery", "Esri.WorldTopoMap")`.
#' @param add_legend Logical to produce a legend or not. Defaults to `TRUE`. Also controls if raster value query appears. 
#' @param add_image_query Logical to add mouse hover query of raster values. Defaults to `TRUE`. 
#' @param max_bytes Maximum size of the raster image in bytes. Defaults to 10MB (10*1024*1024).
#'
#' @return A `leaflet` map object with the data layer(s) and controls.
#'
#' @details
#' - Raster data is projected to WGS84 (EPSG:4326) and colorized using a continuous palette.
#' - Vector data supports POINT, LINESTRING, and POLYGON geometries.
#' - Attributes are displayed in scrollable popups if there are many fields.
#' - The function adds scale bars, measurement tools, and layer controls.
#' - add_image_query can make the file size of the leaflet html widget very large
#'   e.g. up to 20x's larger. It is recommended to have it = `FALSE` when exporting 
#'   as a Web Page or html widget.
#'
#' @examples
#' \dontrun{
#' library(leaflet)
#' library(sf)
#' library(raster)
#' # From scratch
#' leafmap(data = st_read(system.file("shape/nc.shp", package = "sf")))
#'
#' # Add to existing map
#' m <- leaflet()
#' leafmap(m, st_read(system.file("shape/nc.shp", package = "sf")))
#' }
#'
#' @import leaflet sf raster terra
#' @export

leafmap <- function(m = NULL,
                    data = NULL,
                    group_layers = NULL,
                    label = NULL,
                    opacity = 0.5,
                    color = "black",
                    fill_color = color,
                    radius = 3,
                    weight = 2,
                    palette = "viridis",
                    basemaps = c("Esri.WorldImagery", "Esri.WorldTopoMap"),
                    add_legend = TRUE,
                    add_image_query = TRUE,
                    max_bytes = 10 * 1024 * 1024) { # Increased default limit
  
  # —— AUTO‐SWAP FOR STANDALONE VS PIPE —— #
  if (is.null(data) && inherits(m, c("Raster", "SpatRaster", "sf"))) {
    if (is.null(label)) label = paste(substitute(m))
    data <- m
    m <- NULL
  }
  
  stopifnot(inherits(data, c("Raster", "SpatRaster", "sf")))
  
  if (is.null(label)) label = paste(substitute(data))
  
  # Initialize map if needed
  if (is.null(m)) {
    m <- leaflet()
    for (b in basemaps) {
      m <- m %>% addProviderTiles(providers[[b]], group = b)
    }
    group_layers <- NULL
  } else if (is.null(group_layers) && !is.null(attr(m, "group_layers"))) {
    group_layers <- attr(m, "group_layers")
  }
  
  group_layers <- unique(c(group_layers, label))
  
  # —— RASTER PROCESSING —— #
  if (inherits(data, "Raster") || inherits(data, "SpatRaster")) {
    
    # Convert terra to raster for compatibility with projectRasterForLeaflet
    if(inherits(data, "SpatRaster")) data <- raster::raster(data)
    
    # OPTIONAL: Downsample if extremely large to prevent browser crash
    # 2 million cells is a safe threshold for most modern browsers
    if (raster::ncell(data) > 2e6) {
      warning("Raster is very large. Downsampling for performance.")
      data <- raster::aggregate(data, fact = 2, fun = mean)
    }
    
    # Project for Leaflet
    if(is.list(palette) && all(c("classes", "colors") %in% names(palette))){
      sim_leaflet <- round(projectRasterForLeaflet(data, method = "ngb"), 3)
    } else {
      sim_leaflet <- round(projectRasterForLeaflet(data, method = "bilinear"), 3)
    }
    
    raster_vals <- raster::getValues(sim_leaflet)
    
    # Setup Palette
    if (is.list(palette) && all(c("classes", "colors") %in% names(palette))) {
      # Categorical
      pal_classes <- palette$classes
      pal_colors <- palette$colors
      pal_labels <- if (!is.null(palette$labels)) palette$labels else as.character(pal_classes)
      
      pal <- colorFactor(palette = pal_colors, domain = pal_classes, na.color = "#FF000000")
      
      m <- m %>%
        leaflet::addRasterImage(sim_leaflet, colors = pal, opacity = opacity,
                                project = FALSE, layerId = label, group = label,
                                maxBytes = max_bytes) # Fix applied here
    } else {
      # Continuous
      raster_range <- c(raster::minValue(sim_leaflet), raster::maxValue(sim_leaflet))
      if (raster_range[1] < 0 && raster_range[2] > 0) {
        max_abs <- max(abs(raster_range))
        raster_range <- c(-max_abs, max_abs)
      }
      
      pal <- colorNumeric(palette, domain = raster_range, na.color = "#FF000000")
      
      m <- m %>%
        leaflet::addRasterImage(sim_leaflet, colors = pal, opacity = opacity,
                                project = FALSE, layerId = label, group = label,
                                maxBytes = max_bytes) # Fix applied here
    }
    
    # Legend and Query
    if(add_legend){
      if (is.list(palette)){
        m <- leaflet::addLegend(m, colors = pal_colors, labels = pal_labels, title = label, group = label)
      } else {
        m <- leaflet::addLegend(m, pal = pal, values = raster_vals, title = label, group = label,
                                labFormat = labelFormat(big.mark = ""))
      }
    }
    
    if(add_image_query){
      m <- m %>% leafem::addImageQuery(sim_leaflet, project = TRUE, layerId = label, prefix = "")
    }
    
    # JS for Legend Toggling
    m <- m %>%
      htmlwidgets::onRender(sprintf("
        function(el, x) {
          var legend = document.querySelectorAll('.leaflet-control .leaflet-control-legend')[0];
          if (legend) legend.style.display = 'none';
          var map = this;
          map.on('overlayadd', function(e) {
            if (e.name === '%s') { if (legend) legend.style.display = 'block'; }
          });
          map.on('overlayremove', function(e) {
            if (e.name === '%s') { if (legend) legend.style.display = 'none'; }
          });
        }
      ", label, label))
    
    # —— VECTOR PROCESSING (SF) —— #
  } else if (inherits(data, "sf")) {
    x_longlat <- sf::st_transform(data, '+proj=longlat +datum=WGS84')
    geom_type <- unique(sf::st_geometry_type(x_longlat))
    attrs <- sf::st_drop_geometry(data)
    
    popup_content <- sapply(seq_len(nrow(attrs)), function(i) {
      row <- as.list(attrs[i, , drop = FALSE])
      scroll_div_style <- if (ncol(attrs) > 8) "scroll-box" else NULL
      style_block <- if (!is.null(scroll_div_style)) 
        "<style>.scroll-box { max-height:350px; overflow-y:auto; }</style>" else ""
      
      paste0(style_block,
             "<div style='text-align:center;'><strong>", label, "</strong></div><br>",
             "<div", if (!is.null(scroll_div_style)) paste0(" class='", scroll_div_style, "'"), ">",
             "<table style='width:100%;'>",
             paste0("<tr><td><strong>", names(row), "</strong></td><td>", row, "</td></tr>", collapse = ""),
             "</table></div>")
    })
    
    if (any(geom_type %in% c("POINT", "MULTIPOINT"))) {
      m <- m %>% addCircleMarkers(data = x_longlat, radius = radius, weight = weight, color = color, 
                                  fillColor = fill_color, fillOpacity = opacity, popup = popup_content, group = label)
    } else if (any(geom_type %in% c("LINESTRING", "MULTILINESTRING"))) {
      m <- m %>% addPolylines(data = x_longlat, weight = weight, color = color, opacity = opacity, popup = popup_content, group = label)
    } else if (any(geom_type %in% c("POLYGON", "MULTIPOLYGON"))) {
      m <- m %>% addPolygons(data = x_longlat, weight = weight, color = color, fillColor = fill_color, fillOpacity = opacity, popup = popup_content, group = label)
    }
  }
  
  # Final Controls
  m <- m %>%
    leaflet::addLayersControl(baseGroups = basemaps, overlayGroups = group_layers, position = "topleft") %>%
    leaflet::addScaleBar("bottomleft") %>%
    leaflet::addMeasure("bottomleft", primaryLengthUnit = "meters", primaryAreaUnit = "hectares")
  
  attr(m, "group_layers") <- group_layers
  return(m)
}