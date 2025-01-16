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
densityRasterIntersection = function(densA,densB,verbose=F) {
  ## this function locates intersections between rasters
  library(sp)
  library(raster)
  densA <- raster(nrows=40, ncols=40, xmn=0, xmx=2, ymn=0, ymx=2)
  densA[] <- seq(1, 100, length.out=ncell(densA))
  densB <- raster(outer(1:20,20:1), xmn=0, xmx=1, ymn=0, ymx=1)
  densB <- reclassify(densB, c(0,100,NA))
  plot(densA)
  overlapRas = mask(crop(densB, densA), densA)
  plot(overlapRas)
  overlapRas = mask(crop(densA, densB), densB)
  plot(overlapRas)

  return(polygon)
}
