#' @import raster
#' @import MASS
#' @import spocc
#' @import dplyr
#' @import sp
#' @import viridis
#' @import caret
#' @import sf
#' @import rebird
NULL
#' Convert Density Maps to Polygons
#'
#' This function converts density maps (or other raster files) to a polygon
#'
#' @param densityMap Raster of densities from subspeciesDensityMap() or similar
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs_sin = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quant=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' densPol_sin = densityMapToPolygons(densityMap=dens_sin)
densityMapToPolygons = function(densityMap,verbose=F) {
  ## this function converts density maps to polygons
  ## will work on other kinds of polygons as well
  #library(sp)
  #library(raster)
  polygonR = densityMap
  polygonR[!is.na(polygonR)] = 1
  if(verbose==T){print(class(polygonR))}
  polygon = raster::rasterToPolygons(polygonR,fun = NULL,na.rm = T,dissolve = T)
  if(verbose==T){print(class(polygon))}
  polygon = terra::vect(polygon)
  #polygon <- sp::disaggregate(polygon)
  #polygon = st_cast(st_cast,"POLYGON")
  polygon = terra::disagg(polygon)
  polygon = as(sf::st_as_sf(polygon),"Spatial")
  if(verbose==T){print("DISAGGREGATE SUCCESSFUL")}
  #plot(polygon)
  return(polygon)
}
