# load packages
library(runoutSim)
library(terra)
library(sf)

# Load digital elevation model (DEM)
dem <- rast("C:/GitProjects/runoutSIM/Dev/Data/elev_fillsinks_WangLiu.tif")

# Compute hillshade for visualization 
slope <- terrain(dem, "slope", unit="radians")
aspect <- terrain(dem, "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270)

# Load debris flow runout source points and polygons
source_points <- st_read("C:/GitProjects/runoutSIM/Dev/Data/debris_flow_source_points.shp")
runout_polygons <- st_read("C:/GitProjects/runoutSIM/Dev/Data/debris_flow_runout_polygons.shp")

# Load basin boundary
bnd_catchment <- st_read("C:/GitProjects/runoutSIM/Dev/Data/basin_rio_olivares.shp")


map <- leafmap(bnd_catchment, color = '#f7f9f9', fill_color = '#FF000000', weight = 4) %>%
  leafmap(runout_polygons) %>%
  leafmap(source_points, color = '#e74c3c') %>%
  leafmap(hill, palette = grey(0:100/100), opacity = 1, add_legend = FALSE, add_image_query = FALSE) %>%
  leafmap(dem, palette = viridis::mako(10), opacity = 0.6, add_image_query = FALSE)
map

library(leaflet)
library(leafem)

runout_lf <- st_transform(runout_polygons, '+proj=longlat +datum=WGS84')
source_lf <- st_transform(source_points, '+proj=longlat +datum=WGS84')

dem_lf <- round(raster::raster(projectRasterForLeaflet(dem, method = "bilinear" )))

m <- leaflet() %>% addProviderTiles('Esri.WorldImagery') %>%
  #addRasterImage(hill, colors = grey(0:100/100)) %>%
  leaflet::addRasterImage(dem_lf, colors = viridis::viridis(10), opacity = 0.5, layerId = "DEM", group = "DEM", project = FALSE) %>%
  leafem::addImageQuery(dem_lf, layerId = "DEM", project = T) %>%
  addPolygons(data = runout_lf, weight = 1, color = "black", fillColor = "black", fillOpacity = 0.5) %>%
  addCircleMarkers(data = source_lf)


# !! It is leaflem::addImageQuery that is creating the high file size. From 4 MB to 30 to 60 MB depending
# on the rounding of the projectRasterForLeaflet raster...
