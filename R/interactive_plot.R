require(leaflet)
require(leafem)
require(raster)


# ! Need to update to make flexible to multiple rasters and vectors
# Rasters = source probability, runout walks, connectivity probability, source cells
# Vector = source points, source polygons, runout polygons, connectivity feature

# Plot requirements
# a. source points (vector)
# b. runout polygons (vector)
# c. connect feature polygon (vector)
# d. connectivity (connect probability) (raster)
# e. runout simulation (traverse probability) (raster)

plot.runout <- function(x, y = NULL, z = NULL, h = NULL, 
                        rasterTitle = "Traverse Probability",
                        polyTitle = "Observed Runout", 
                        rasterOpacity = 0.9,
                        polyOpacity = 0.1,
                        rasterPalette = "viridis",
                        polyColor = "black"){
  
  # Create color palette
  pal <- colorNumeric(rasterPalette, 0:1, na.color = "#FF000000")
  
  # Convert the raster to a leaflet-compatible format
  sim_leaflet <- round(raster::raster(projectRasterForLeaflet(x, method = "bilinear")),3)
  
  sim.leaflet_map <- leaflet() %>%
    
    addProviderTiles(providers$Esri.WorldImagery)
  
  group_layers = NULL
  
  if(!is.null(y)){
    
    poly_leaflet <- st_transform(st_as_sf(y), '+proj=longlat +datum=WGS84')
    
    sim.leaflet_map <- sim.leaflet_map %>%
      
      addPolygons(data = poly_leaflet, group = polyTitle,
                  color = polyColor, fillOpacity = polyOpacity,
                  weight = 3, opacity = 0.9)
    
    group_layers = polyTitle
    
    
  }
  
  
  if(!is.null(z)){
    
    # for map of connectivity
    z_leaflet <- round(raster::raster(projectRasterForLeaflet(z, method = "bilinear")),3)
    
    sim.leaflet_map <- sim.leaflet_map %>%
      addRasterImage(z_leaflet, colors = pal, opacity = rasterOpacity, 
                     project = TRUE, group = "Connect Probability", layerId = "Connect Probability") %>%
      addImageQuery(z_leaflet, project = TRUE, layerId  = "Connect Probability", prefix = "") %>%
      addLegend(pal = pal, values=0:1, title = "Connect Probability")
    
    group_layers <- c(group_layers, "Connect Probability")
    
  }
  
  # for map of runout
  sim.leaflet_map <- sim.leaflet_map %>% 
    addRasterImage(sim_leaflet, colors = pal, opacity = rasterOpacity, 
                   project = TRUE, layerId  = rasterTitle, group = rasterTitle) %>%
    addImageQuery(sim_leaflet, project = TRUE, layerId  = rasterTitle, prefix = "") # project = TRUE for this to work.
  
  group_layers <- c(group_layers, rasterTitle)
  
  sim.leaflet_map <- sim.leaflet_map %>%    
    
    addLegend(pal = pal, values=0:1, title = rasterTitle) %>%
    
    addLayersControl(overlayGroups = group_layers,
                     options = layersControlOptions(collapsed = FALSE)) %>%
    
    addScaleBar("bottomleft") %>%
    
    addMeasure("bottomleft", 
               primaryLengthUnit = "meters", 
               secondaryLengthUnit = "kilometers",
               primaryAreaUnit = "hectares",
               secondaryAreaUnit = "sqmeters")
  
  return(sim.leaflet_map)
  
}

