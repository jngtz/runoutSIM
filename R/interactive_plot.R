#' Plot Spatial Data Using Leaflet
#'
#' Creates an interactive `leaflet` map from raster or `sf` vector data, with options to customize basemaps, styling, and popups.
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
#' @param palette Color palette name used with `leaflet::colorNumeric()` for raster coloring. Defaults to `"viridis"`. If categorical values, supply a list - e.g. list(classes = 1, colors = "#99d2ff", labels = "Water"))
#' @param basemaps Character vector of tile provider names (from `leaflet::providers`) to include as base layers. Defaults to `c("Esri.WorldImagery", "Esri.WorldTopoMap")`.
#' @param add_legend Logical to produce a legend or not. Defaults to `TRUE`.
#' @param add_image_query Logical to add mouse hover query of raster values. Defaults to `TRUE`. Independent of `add_legend`.
#'
#' @return A `leaflet` map object with the data layer(s) and controls.
#'
#' @import leaflet sf raster terra leafem htmlwidgets
#' @importFrom magrittr %>%
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
                    add_image_query = TRUE) {
  
  # —— AUTO-SWAP FOR STANDALONE VS PIPE —— #
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
    m <- leaflet::leaflet()
    for (b in basemaps) {
      m <- m %>% leaflet::addProviderTiles(leaflet::providers[[b]], group = b)
    }
    group_layers <- NULL
  } else if (is.null(group_layers) && !is.null(attr(m, "group_layers"))) {
    group_layers <- attr(m, "group_layers")
  }
  
  # Helper to track groups
  group_layers <- unique(c(group_layers, label))
  
  # If raster
  if (inherits(data, "Raster") || inherits(data, "SpatRaster")) {
    
    if (inherits(data, "SpatRaster")){
      data <- raster::raster(data)
    }
    
    if (is.list(palette) && all(c("classes", "colors") %in% names(palette))){
      sim_leaflet <- round(leaflet::projectRasterForLeaflet(data, method = "ngb"), 3)
    } else {
      sim_leaflet <- round(leaflet::projectRasterForLeaflet(data, method = "bilinear"), 3)
    }
    
    raster_vals <- raster::getValues(sim_leaflet)
    
    # Unique class so this layer's legend can be targeted individually in JS below
    # (fixes: previously all legends were selected via a hardcoded [0] index)
    legend_class <- paste0("legend-", gsub("[^A-Za-z0-9]", "", label))
    
    if (is.list(palette) && all(c("classes", "colors") %in% names(palette))) {
      # Handle categorical raster
      pal_classes <- palette$classes
      pal_colors <- palette$colors
      pal_labels <- if (!is.null(palette$labels)) palette$labels else as.character(pal_classes)
      
      if (length(pal_classes) != length(pal_colors)) stop("`palette$classes` and `palette$colors` must have the same length.")
      if (length(pal_labels) != length(pal_classes)) stop("`palette$labels` must match the number of classes.")
      
      pal <- leaflet::colorFactor(palette = pal_colors, domain = pal_classes, na.color = "#FF000000")
      
      m <- m %>%
        leaflet::addRasterImage(sim_leaflet, colors = pal, opacity = opacity,
                                project = FALSE, layerId = label, group = label) 
      
      if (add_image_query){
        m <- m %>% 
          leafem::addImageQuery(sim_leaflet, project = TRUE, layerId = label, prefix = "")
      }
      
      if (add_legend){
        m <- leaflet::addLegend(m,
                                colors = pal_colors,
                                labels = pal_labels,
                                title = label,
                                group = label,
                                className = paste("info legend", legend_class)
        )
      }
      
    } else {
      # Handle continuous raster
      raster_range <- c(raster::minValue(sim_leaflet), raster::maxValue(sim_leaflet))
      
      if (raster_range[1] < 0 && raster_range[2] > 0) {
        max_abs <- max(abs(raster_range))
        raster_range <- c(-max_abs, max_abs)
      }
      
      pal <- leaflet::colorNumeric(palette, domain = raster_range, na.color = "#FF000000")
      
      m <- m %>%
        leaflet::addRasterImage(sim_leaflet, colors = pal, opacity = opacity,
                                project = FALSE, layerId = label, group = label) 
      
      if (add_legend){
        m <- leaflet::addLegend(m, pal = pal, values = raster_vals, title = label, group = label,
                                labFormat = leaflet::labelFormat(big.mark = ""),
                                className = paste("info legend", legend_class))
      }
      
      # add_image_query is independent of add_legend now (bug fix: previously
      # nested inside `if (add_legend)`, so turning off the legend silently
      # disabled hover-query too, unlike the categorical branch above)
      if (add_image_query){
        m <- m %>% 
          leafem::addImageQuery(sim_leaflet, project = TRUE, layerId = label, prefix = "")
      }
    }
    
    # Show/hide legend with layer toggles
    # (bug fix: target this layer's legend specifically via legend_class,
    # instead of always grabbing index [0] which broke with multiple raster layers)
    m <- m %>%
      htmlwidgets::onRender(sprintf("
        function(el, x) {
          var legend = document.querySelector('.leaflet-control.%s');
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
      ", legend_class, label, label))
    
  } else if (inherits(data, "sf")) {
    x_longlat <- sf::st_transform(data, '+proj=longlat +datum=WGS84')
    geom_type <- unique(sf::st_geometry_type(x_longlat))
    
    attrs <- sf::st_drop_geometry(data)
    
    popup_content <- sapply(seq_len(nrow(attrs)), function(i) {
      row <- as.list(attrs[i, , drop = FALSE])
      n_rows <- nrow(attrs)
      scroll_div_style <- if (n_rows > 8) "scroll-box" else NULL
      
      style_block <- if (!is.null(scroll_div_style)) paste0(
        "<style>",
        ".", scroll_div_style, " { max-height:350px; overflow-y:auto; }",
        ".", scroll_div_style, "::-webkit-scrollbar { width:6px; }",
        ".", scroll_div_style, "::-webkit-scrollbar-track { background:transparent; border-radius:10px; }",
        ".", scroll_div_style, "::-webkit-scrollbar-thumb { background:rgba(0,0,0,0.3); border-radius:10px; }",
        ".", scroll_div_style, " { scrollbar-width: 5px; scrollbar-color: rgba(0,0,0,0.3) transparent; }",
        "</style>"
      ) else ""
      
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
    
    if (any(geom_type %in% c("POINT", "MULTIPOINT"))) {
      m <- m %>% leaflet::addCircleMarkers(
        data = x_longlat,
        radius = radius,
        weight = weight,
        color = color,
        fillColor = fill_color,
        fillOpacity = opacity,
        # bug fix: `~label` was a formula referring to a data column named
        # "label" (which almost never exists); use the `label` string itself
        label = label,
        popup = popup_content,
        group = label
      )
    } else if (any(geom_type %in% c("LINESTRING", "MULTILINESTRING"))) {
      m <- m %>% leaflet::addPolylines(
        data = x_longlat,
        weight = weight,
        color = color,
        opacity = opacity,
        popup = popup_content,
        group = label
      )
    } else if (any(geom_type %in% c("POLYGON", "MULTIPOLYGON"))) {
      m <- m %>% leaflet::addPolygons(
        data = x_longlat,
        weight = weight,
        color = color,
        fillColor = fill_color,
        fillOpacity = opacity,
        popup = popup_content,
        group = label
      )
    } else {
      warning(paste("Unsupported geometry type:", paste(geom_type, collapse = ", ")))
    }
  } else {
    warning("Input object 'data' must be a raster or sf object.")
  }
  
  # Update map with controls and tools
  m <- m %>%
    leaflet::addLayersControl(
      baseGroups = basemaps,
      overlayGroups = group_layers,
      position = "topleft",
      options = leaflet::layersControlOptions(collapsed = TRUE)
    ) %>%
    leaflet::addScaleBar("bottomleft") %>%
    leaflet::addMeasure("bottomleft",
                        primaryLengthUnit = "meters",
                        secondaryLengthUnit = "kilometers",
                        primaryAreaUnit = "hectares",
                        secondaryAreaUnit = "sqmeters")
  
  attr(m, "group_layers") <- group_layers
  return(m)
}