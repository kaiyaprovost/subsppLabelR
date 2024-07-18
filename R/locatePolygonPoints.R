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
#' Point to Polygon Locator
#'
#' This function takes a one polygon and a list of points
#' and determines whether the points fall within the polygon, then
#' adds a column for each polygon's name with a boolean of whether
#' the point overlaps.
#'
#' @param test_points The points to check where they fall
#' @param polygonA The  polygon
#' @param setcoord Whether to set a coordinate reference system to the points and polygons
#' @param crs A coordinate reference system to assign to the polygons and the test points if setcoord = TRUE
#' @param nameA The (subspecies) name associated with the first polygon
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_sin = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_ful = labeledLoc[labeledLoc$subspecies=="fulvescens",]
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quant=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' dens_ful = subspeciesDensityMap(localities=locs_ful,quant=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' densPol_sin = densityMapToPolygons(densityMap=dens_sin)
#' densPol_ful = densityMapToPolygons(densityMap=dens_ful)
#' polyLocations = labeledLoc
#' polyLocations = locatePolygonPoints(test_points=polyLocations,polygonA=densPol_sin,
#'    crs="+proj=longlat +ellps=WGS84",setcoord = TRUE,nameA="sinuatus")
locatePolygonPoints = function(test_points,
                                polygonA,
                                crs = "+proj=longlat +ellps=WGS84",
                                setcoord = TRUE,
                                nameA = "A") {
  pts = test_points
  pts$longitude = as.numeric(pts$longitude)
  pts$latitude = as.numeric(pts$latitude)
  if (setcoord == T) {
    sp::coordinates(pts) = ~ longitude + latitude
    sp::proj4string(pts) = sp::CRS(crs)
    polygonA = sp::spTransform(polygonA, crs)
  }
  overpolygon = cbind(sp::over(pts, polygonA))
  overpolygon[is.na(overpolygon)] = 0
  colnames(overpolygon) = nameA
  test_points = cbind(test_points,overpolygon)
  return(test_points)
}
