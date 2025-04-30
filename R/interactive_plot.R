require(leaflet)
require(leafem)

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
#'
#' @return A `leaflet` map object with the data layer(s) and controls.
#'
#' @details
#' - Raster data is projected to WGS84 (EPSG:4326) and colorized using a continuous palette.
#' - Vector data supports POINT, LINESTRING, and POLYGON geometries.
#' - Attributes are displayed in scrollable popups if there are many fields.
#' - The function adds scale bars, measurement tools, and layer controls.
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

leafmap<- function(m = NULL,
                        data = NULL,               # ← default data to NULL
                        group_layers = NULL,
                        label = NULL,
                        opacity = 0.5,
                        color = "black",
                        fill_color = color,
                        radius = 3,
                        weight = 2,
                        palette = "viridis",
                        basemaps = c("Esri.WorldImagery", "Esri.WorldTopoMap"),
                        add_legend = TRUE) {
  
  # —— AUTO‐SWAP FOR STANDALONE VS PIPE —— #
  # if data was never given, but m *is* a raster or sf, assume user called
  # plot.runout(raster, …) and swap them
  if (is.null(data) && inherits(m, c("Raster", "SpatRaster", "sf"))) {
    if (is.null(label)){
      label = paste(substitute(m))
    }
    data <- m
    m <- NULL
  }
  
  stopifnot(inherits(data, c("Raster", "SpatRaster", "sf")))
  
  if (is.null(label)){
    label = paste(substitute(data))
  }
  
  # Initialize map if needed
  if (is.null(m)) {
   m <- leaflet()
    # add each provider as a base layer
    for (b in basemaps) {
      m <- m %>% addProviderTiles(providers[[b]], group = b)
    }
    group_layers <- NULL
  } else if (is.null(group_layers) && !is.null(attr(m, "group_layers"))) {
    group_layers <- attr(m, "group_layers")
  }
  
  # Helper to track groups
  group_layers <- unique(c(group_layers, label))
  
  # if raster
  if (inherits(data, "Raster") || inherits(data, "SpatRaster")) {
    
    if(inherits(data, "SpatRaster")){
      data <- raster::raster(data)
    }
    
    if(is.list(palette) && all(c("classes", "colors") %in% names(palette))){
      sim_leaflet <- round(projectRasterForLeaflet(data, method = "ngb"), 3)
    } else {
      sim_leaflet <- round(projectRasterForLeaflet(data, method = "bilinear"), 3)
    }
    
    
    raster_vals <- raster::getValues(sim_leaflet)
    
    if (is.list(palette) && all(c("classes", "colors") %in% names(palette))) {
      # Handle categorical raster
      pal_classes <- palette$classes
      pal_colors <- palette$colors
      pal_labels <- if (!is.null(palette$labels)) palette$labels else as.character(pal_classes)
      
      if (length(pal_classes) != length(pal_colors)) stop("`palette$classes` and `palette$colors` must have the same length.")
      if (length(pal_labels) != length(pal_classes)) stop("`palette$labels` must match the number of classes.")
      
      pal <- colorFactor(palette = pal_colors, domain = pal_classes, na.color = "#FF000000")
      
      m <- m %>%
        leaflet::addRasterImage(sim_leaflet, colors = pal, opacity = opacity,
                       project = TRUE, layerId = label, group = label) %>%
        leafem::addImageQuery(sim_leaflet, project = TRUE, layerId = label, prefix = "")
      
      if(add_legend){
        m <- leaflet::addLegend(m,
                                colors = pal_colors,
                                labels = pal_labels,
                                title = label,
                                group = label
        )
      }
      

      
    } else {
      # Handle continuous raster
      #raster_range <- range(raster_vals, na.rm = TRUE)
      raster_range <- c(raster::minValue(sim_leaflet), raster::maxValue(sim_leaflet))
      
      # Ensure symmetric range around zero if values span negative and positive
      if (raster_range[1] < 0 && raster_range[2] > 0) {
        max_abs <- max(abs(raster_range))
        raster_range <- c(-max_abs, max_abs)
      }
      
      
      pal <- colorNumeric(palette, domain = raster_range, na.color = "#FF000000")
      
      m <- m %>%
        leaflet::addRasterImage(sim_leaflet, colors = pal, opacity = opacity,
                       project = TRUE, layerId = label, group = label) 
        
        
      if(add_legend){
        m <- leaflet::addLegend(m, pal = pal, values = raster_vals, title = label, group = label,
                           labFormat = labelFormat(big.mark = "")) %>%
          leafem::addImageQuery(sim_leaflet, project = TRUE, layerId = label, prefix = "")
      }

    }
    
    # Show/hide legend with layer toggles
    m <- m %>%
      htmlwidgets::onRender(sprintf("
        function(el, x) {
          var legend = document.querySelectorAll('.leaflet-control .leaflet-control-legend')[0];
          if (legend) legend.style.display = 'none';

          var map = this;
          map.on('overlayadd', function(e) {
            if (e.name === '%s') {
              if (legend) legend.style.display = 'block';
            }
          });
          map.on('overlayremove', function(e) {
            if (e.name === '%s') {
              if (legend) legend.style.display = 'none';
            }
          });
        }
      ", label, label))
    
  } else if (inherits(data, "sf")) {
    x_longlat <- sf::st_transform(data, '+proj=longlat +datum=WGS84')
    geom_type <- unique(sf::st_geometry_type(x_longlat))
    
    # Create popup content from attributes
    attrs <- sf::st_drop_geometry(data)
    
    popup_content <- sapply(seq_len(nrow(attrs)), function(i) {
      row     <- as.list(attrs[i, , drop = FALSE])
      n_rows  <- nrow(attrs)
      scroll_div_style <- if (n_rows > 8) {
        # class name on the wrapper; actual CSS lives in the <style> block
        "scroll-box"
      } else {
        NULL
      }
      
      # 1) Define a little <style> block up front
      style_block <- if (!is.null(scroll_div_style)) paste0(
        "<style>",
          ".", scroll_div_style, " { max-height:350px; overflow-y:auto; }",
          ".", scroll_div_style, "::-webkit-scrollbar { width:6px; }",
        ".", scroll_div_style, "::-webkit-scrollbar-track { background:transparent; border-radius:10px; }",
        ".", scroll_div_style, "::-webkit-scrollbar-thumb { background:rgba(0,0,0,0.3); border-radius:10px; }",
          ".", scroll_div_style, " { scrollbar-width: 5px; scrollbar-color: rgba(0,0,0,0.3) transparent; }",
        "</style>"
      ) else ""
      
      # 2) Wrap the table in a <div class='scroll-box'>…</div>
      paste0(
        style_block,
        "<div style='text-align:center;'><strong style='font-size:14px;'>", label, "</strong></div><br>",
        "<div", if (!is.null(scroll_div_style)) paste0(" class='", scroll_div_style, "'"), ">",
        "<table style='border-collapse:collapse; width:100%; text-align:left;'>",
        paste0(
          "<tr><td style='padding:2px 4px;'><strong>", names(row), "</strong></td>",
          "<td style='padding:2px 4px;'>", row, "</td></tr>",
          collapse = ""
        ),
        "</table>",
        "</div>"
      )
    })
    
    if (all(geom_type %in% c("POLYGON", "MULTIPOLYGON"))) {
      m <- m %>%
        leaflet::addPolygons(data = x_longlat, group = label,
                    color = color, fillColor = fill_color, fillOpacity = opacity,
                    weight = weight, opacity = 0.9,
                    popup = popup_content[1])
      
    } else if (all(geom_type %in% c("LINESTRING", "MULTILINESTRING"))) {
      m <- m %>%
        leaflet::addPolylines(data = x_longlat, group = label,
                     color = color, weight = weight, opacity = opacity,
                     popup = popup_content)
      
    } else if (all(geom_type %in% c("POINT", "MULTIPOINT"))) {
      m <- m %>%
        leaflet::addCircleMarkers(data = x_longlat, group = label,
                         radius = radius, color = color, fillColor = fill_color,
                         stroke = TRUE, fillOpacity = opacity,
                         popup = ~popup_content)
    } else {
      warning(paste("Unsupported geometry type:", paste(geom_type, collapse = ", ")))
    }
  } else {
    warning("Input object 'data' must be a raster or sf object.")
  }
  
  # Update map with controls and tools
  m <- m %>%
    leaflet::addLayersControl(
      baseGroups    = basemaps,
      overlayGroups = group_layers,
      position = "topleft",
      options       = layersControlOptions(collapsed = TRUE)
    ) %>%
    leaflet::addScaleBar("bottomleft") %>%
    leaflet::addMeasure("bottomleft",
               primaryLengthUnit   = "meters",
               secondaryLengthUnit = "kilometers",
               primaryAreaUnit     = "hectares",
               secondaryAreaUnit   = "sqmeters")
  
  attr(m, "group_layers") <- group_layers
  return(m)

}

