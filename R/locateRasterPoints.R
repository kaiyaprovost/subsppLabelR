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
#' Point to Raster Locator
#'
#' This function takes one raster and a list of points
#' and determines whether the points fall within the raster, then
#' adds a column for each raster's name with a boolean of whether
#' the point overlaps.
#'
#' @param test_points The points to check where they fall
#' @param rasterA The Raster to look for overlap  
#' @param nameA The (subspecies) name associated with the raster
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
#' polyLocations = locateRasterPoints(test_points=polyLocations,polygonA=densPol_sin,
#'    nameA="sinuatus")
locateRasterPoints = function(test_points,
                                rasterA,
                                nameA = "A") {
  pts = test_points
  pts$longitude = as.numeric(pts$longitude)
  pts$latitude = as.numeric(pts$latitude)
  existingCells=cellFromXY(rasterA,pts[,c("longitude","latitude")])
  valuesCells = values(rasterA)[existingCells]
  valuesCells[!(is.na(valuesCells))] = 1
  valuesCells[(is.na(valuesCells))] = 0
  valuesCells = as.data.frame(valuesCells)
  colnames(valuesCells) = nameA
  test_points = cbind(test_points,valuesCells)
  return(test_points)
}
