require(leaflet)
require(leafem)


plot.leaflet <- function(m = NULL,
                        data = NULL,               # ← default data to NULL
                        group_layers = NULL,
                        label = NULL,
                        opacity = 0.5,
                        color = "black",
                        radius = 3,
                        weight = 2,
                        palette = "viridis",
                        basemaps = c("Esri.WorldImagery", "Esri.WorldTopoMap")) {
  
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
  
  # If raster
  if (inherits(data, "Raster") || inherits(data, "SpatRaster")) {
    # Project and round raster for display
    sim_leaflet <- round(raster::raster(projectRasterForLeaflet(data, method = "bilinear")), 3)
    
    # Get min/max values from raster
    raster_range <- range(values(sim_leaflet), na.rm = TRUE)
    
    # Create color palette based on actual raster values
    pal <- colorNumeric(palette, domain = raster_range, na.color = "#FF000000")
    
    m <- m %>%
      addRasterImage(sim_leaflet, colors = pal, opacity = opacity,
                     project = TRUE, layerId = label, group = label) %>%
      addImageQuery(sim_leaflet, project = TRUE, layerId = label, prefix = "") %>%
      addLegend(pal = pal, values = values(sim_leaflet), title = label)
  } else if (inherits(data, "sf")) {
    x_longlat <- st_transform(data, '+proj=longlat +datum=WGS84')
    geom_type <- unique(st_geometry_type(x_longlat))
    
    # Create popup content from attributes
    attrs <- st_drop_geometry(data)
    
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
        addPolygons(data = x_longlat, group = label,
                    color = color, fillOpacity = opacity,
                    weight = weight, opacity = 0.9,
                    popup = popup_content[1])
      
    } else if (all(geom_type %in% c("LINESTRING", "MULTILINESTRING"))) {
      m <- m %>%
        addPolylines(data = x_longlat, group = label,
                     color = color, weight = weight, opacity = opacity,
                     popup = popup_content)
      
    } else if (all(geom_type %in% c("POINT", "MULTIPOINT"))) {
      m <- m %>%
        addCircleMarkers(data = x_longlat, group = label,
                         radius = radius, color = color,
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
    addLayersControl(
      baseGroups    = basemaps,
      overlayGroups = group_layers,
      position = "topleft",
      options       = layersControlOptions(collapsed = TRUE)
    ) %>%
    addScaleBar("bottomleft") %>%
    addMeasure("bottomleft",
               primaryLengthUnit   = "meters",
               secondaryLengthUnit = "kilometers",
               primaryAreaUnit     = "hectares",
               secondaryAreaUnit   = "sqmeters")
  
  attr(m, "group_layers") <- group_layers
  return(m)

}

