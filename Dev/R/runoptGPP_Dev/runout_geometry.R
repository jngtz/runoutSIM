# GEOMETRY CORE FUNCTIONS #######################################################

#' Minimum area bounding rectangles
#'
#' Computes a minimum area bounding rectangle for a set of points
#' @param pts A matrix with x and y coordinates of points
#' @return A list including the bounding box coordinates
#' @examples
#' pts=cbind(runif(60),runif(60))
#' bb=minbb(pts)
#' plot(bb$box,type="l")
#' points(pts)
#' pts2=rotxy(pts,pi/1.6)
#' bb2=minbb(pts2)
#' plot(bb2$box,type="l")
#' points(pts2)

minbb <- function(pts){
  #Minimum area bounding rectangles
  #Code thanks to Barry Rowlingson on R-sig-Geo mailing list
  #https://stat.ethz.ch/pipermail/r-sig-geo/2011-March/011212.html

  # compute convex hull of a set of points
  #(giving indices of points lying on the convex hull in clockwise order)
  ch = chull(pts)
  pts=pts[ch,] #make matrix of convex hull point locations (x,y)
  np = nrow(pts) #count the number of points
  pts=rbind(pts,pts[1,]) #add  first point to end of matrix to make polygon
  minbba = Inf ; bbth = NA; rotmin = NA
  for(i in 1:np){
    th = pi-atan2(pts[i+1,2]-pts[i,2],pts[i+1,1]-pts[i,1]) #edge directions
    prot = rotxy(pts,th)
    bba = diff(range(prot[,1])) * diff(range(prot[,2])) #calc box area

    # if bounding box area is less than min bounding box
    if(bba < minbba){
      xyb=cbind(
        c(min(prot[,1]),max(prot[,1]),max(prot[,1]),min(prot[,1])),
        c(min(prot[,2]),min(prot[,2]),max(prot[,2]),max(prot[,2]))
      )
      xyb=rbind(xyb,xyb[1,])
      xyb= rotxy(xyb,-th)

      rotmin=prot
      minbba = bba
      bbth = th
    }
  }
  return(list(minbba=minbba,theta=th,pts=rotmin,box=xyb))

}


#' Rotate bounding box
#'
#' @param pts A matrix with x and y coordinates of points
#' @param angle Angle to rotate
#' @return A list including the bounding box coordinates
#' @examples
#' pts=cbind(runif(60),runif(60))
#' bb=minbb(pts)
#' plot(bb$box,type="l")
#' points(pts)
#' pts2=rotxy(pts,pi/1.6)
#' bb2=minbb(pts2)
#' plot(bb2$box,type="l")
#' points(pts2)
#' rotxy(pts, 30)

rotxy = function(pts,angle){
  #Rotates bounding box
  co = cos(angle)
  si = sin(angle)
  cbind(co * pts[,1] - si * pts[,2], si * pts[,1] + co * pts[,2])
}



#' Calculate Euclidean distance between two points
#'
#' @param p1x x coordinate of point 1
#' @param p1y y coordinate of point 1
#' @param p2x x coordinate of point 2
#' @param p2y y coordinate of point 2
#' @return The distance between two points
#' @examples
#' EuclDist(2,3,7,8)

EuclDist <- function(p1x, p1y, p2x, p2y){
  #Calculate Euclidean distance between two points
  dist <- sqrt((p2x - p1x)^2 + (p2y - p1y)^2)
  return(dist)
}


#' Get vertices from a polygon
#'
#' @param x A SpatialPolygon object
#' @return A matrix with all verticies

getVertices <- function(x){
  # Returns matrix of all the vertices in a polygon
  # Can us as input to minbb()
  n_polys <- length(x@polygons[[1]]@Polygons)
  m_coords <- x@polygons[[1]]@Polygons[[1]]@coords

  if(n_polys > 1){
    for(i in 2:n_polys){
      m_coords <- rbind(m_coords, x@polygons[[1]]@Polygons[[i]]@coords)
    }
  }

  m_coords

}




# APPLIED FUNCTIONS ####################################################

#' Minimum area bounding box for spatial polgyons
#'
#' Determines min. area bounding box for a single or set of spatial polygons
#' @param x A SpatialPolygonsDataFrame
#' @return A SpatialPolygonsDataFrame with corresponding bounding boxes
#' @examples
#' file_nm <- system.file("extdata/dflow_runout_ply.shp", package="runout.opt")
#' slide_plys <- rgdal::readOGR(file_nm)
#' minbb <- minBBoxSpatialPolygons(slide_plys)
#'
#' sp::plot(slide_plys)
#' sp::plot(minbb, add = TRUE)

minBBoxSpatialPolygons <- function(x) {

  #number of features
  n_features <- length(x@polygons)

  bbList <- list()
  for (i in 1:n_features) {
    ply <- x[i,]
    #Get the vertices coordinates of from SpatialPolygons
    #pnts <- ply@polygons[[1]]@Polygons[[1]]@coords
    pnts <- getVertices(ply)
    #Calculate minimum area rectangle
    bb <- minbb(pnts)
    #Create and add to a list of polygons
    bbply <- sp::Polygon(bb$box)
    bbplys <- sp::Polygons(list(bbply), as.character(i))
    bbList[[i]]<- bbplys
    #Output shapefile of bounding boxes
    bbSp <- sp::SpatialPolygons(bbList, 1:length(bbList))
    bbSpDf <- sp::SpatialPolygonsDataFrame(bbSp,
                                       data.frame(value = 1:length(bbList),
                                                  row.names = 1:length(bbList)))
    sp::proj4string(bbSpDf) <- sp::proj4string(x)
  }
  return(bbSpDf)
}




#' Runout Geometry
#'
#' Computes runnout attributes: width, length, area, max.
#'     min. elevation, and surface area using Spatial Polygons and a DEM
#' @param runout_plys A SpatialPolygonsDataFrame of runout tracks
#' @param elev A RasterLayer DEM
#' @param ID NOT SURE...
#' @return A data frame containing runout geometries
#' @examples
#' # Load elevation model (DEM)
#' dem <- raster::raster(system.file("extdata/elev_12_5m.tif", package="runout.opt"))
#'
#' # Load runout polygons
#' file_nm <- system.file("extdata/dflow_runout_ply.shp", package="runout.opt")
#' slide_plys <- rgdal::readOGR(file_nm)
#'
#' slide_geom <- runoutGeom(slide_plys, dem)
#' slide_geom

runoutGeom <- function(runout_plys, elev, ID = NULL) {
  #Create a dataframe to store the ID, length, width and (landslide) area
  
  if(class(runout_plys)[1] == "sf"){
    runout_plys = sf::as_Spatial(runout_plys)
  }
  
  if(is.null(ID)){
    bbDf <- data.frame(fid = 0:(length(runout_plys)-1), id = 1:length(runout_plys),
                       width = NA, length = NA, area = NA, surfacearea = NA,
                       maxelev = NA, minelev = NA, reachangle = NA)
  } else {
    bbDf <- data.frame(fid = 0:(length(runout_plys)-1), id = NA,
                       width = NA, length = NA, area = NA, surfacearea = NA,
                       maxelev = NA, minelev = NA, reachangle = NA)
  }


  n_features <- length(runout_plys@polygons)

  for (i in 1:n_features){
    runout_ply <- runout_plys[i,]

    if(!is.null(ID)){bbDf[i,]$id <- runout_ply@data[,ID]}
    #Get the vertices coordinates of from SpatialPolygons
    #pnts <- runout_ply@polygons[[1]]@Polygons[[1]]@coords
    pnts <- getVertices(runout_ply)
    #Calculate minimum area rectangle
    bb <- minbb(pnts)

    #Calcluate the difference in elevation between points
    elevExt <- raster::extract(elev, bb$box[1:4,])
    dElev12 <- sqrt( (elevExt[1] - elevExt[2])^2)
    dElev14 <- sqrt( (elevExt[1] - elevExt[4])^2)

    if(is.na(dElev14)){
      #If bounding box node extends beyond dem (e.g. in NA area)
      dElev14 <- sqrt( (elevExt[2] - elevExt[3])^2)
    }
    #Determine (and calculate) the length and width based on delta elevation

    if(is.na(dElev12)){
      #If bounding box node extends beyond dem (e.g. in NA area)
      dElev12 <- sqrt( (elevExt[4] - elevExt[3])^2)
    }

    if(dElev12 > dElev14) {
      length <- EuclDist(bb$box[1,1], bb$box[1,2], bb$box[2,1], bb$box[2,2])
      width <- EuclDist(bb$box[1,1], bb$box[1,2], bb$box[4,1], bb$box[4,2])
    } else {
      length <- EuclDist(bb$box[1,1], bb$box[1,2], bb$box[4,1], bb$box[4,2])
      width <- EuclDist(bb$box[1,1], bb$box[1,2], bb$box[2,1], bb$box[2,2])
    }

    bbDf[i,]$length <- length #planar
    bbDf[i,]$width <- width #planar
    #bbDf[i,]$area <- rgeos::gArea(runout_ply) #area of landslide (*not area of bbox) {Rgeos}
    bbDf[i,]$area <- as.numeric(sf::st_area(sf::st_as_sf(runout_ply)))

    #Calculate the 'true' surface area
    elevCrop <- raster::crop(elev, runout_ply) #crop and mask used to speed up calculation
    elevMask <- raster::mask(elevCrop, runout_ply)
    elevMaskSGDF <- as(elevMask , "SpatialGridDataFrame")
    bbDf[i,]$surfacearea <- sp::surfaceArea(elevMaskSGDF) #surfacearea of landslide {sp}
    #Calculate the max. and min. elevation
    elevPnts <- raster::extract(elev, pnts)
    bbDf[i,]$maxelev <- max(elevPnts)
    bbDf[i,]$minelev <- min(elevPnts)
    bbDf[i,]$reachangle <- 90 - (atan( bbDf[i,]$length / ( bbDf[i,]$maxelev- bbDf[i,]$minelev))*180/pi)

  }
  return(bbDf)
}


