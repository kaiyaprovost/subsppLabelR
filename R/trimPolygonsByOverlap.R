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
#' Remove Overlapping Polygons
#'
#' This function takes polygon features flagged for removal from flagPolygonOverlap()
#' and removes them from the polygon. Called within polygonTrimmer().
#'
#' @param polygon Polygon to remove flagged features from.
#' @param idList List of flagged features to completely remove.
#' @param intList List of intersections to remove.
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
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
#' polygonsToRemove = (flagPolygonOverlap(subsppPoly1=densPol_sin,
#'    subsppPoly2=densPol_ful))
#' sin_polygonsToRemove = polygonsToRemove$subsppApoly_toremove
#' ful_polygonsToRemove = polygonsToRemove$subsppBpoly_toremove
#' overlapToRemove_sin = polygonsToRemove$subsppA_intToRemove
#' overlapToRemove_ful = polygonsToRemove$subsppB_intToRemove
#' sin_densityPolygon_trim = trimPolygonsByOverlap(polygon=densPol_sin,
#'    idList = sin_polygonsToRemove,intList=overlapToRemove_sin)
#' ful_densityPolygon_trim = trimPolygonsByOverlap(polygon=densPol_ful,
#'    idList = ful_polygonsToRemove,intList=overlapToRemove_ful)
trimPolygonsByOverlap = function(polygon,
                                 idList = NULL,
                                 intList = NULL) {
  ## function(polygon,idList=NULL,intList=NULL)
  ## it then removes polygons that are completely within another polygon as ID'd by flagPolygonOverlap()
  ## TODO: what about things that are in neither polygon?
  ## TODO: what about things that are in both?
  polytrim = polygon
  if (is.null(idList)) {
    if (is.null(intList)) {
      #print("no overlaps to remove")
      polytrim = rgeos::gUnion(polytrim, polytrim)
    } else {
      for (int in 1:length(intList)) {
        polytrim = rgeos::gDifference(polytrim, intList[[int]])
      }
      #print("removing differences from larger shapefile")
    }
  } else {
    ## add functionality in case removes the whole polygon?
    if(length(idList)<length(polygon)){
      polytrim = polygon[-idList,]
    } else {
      ## whole polygon is going to be removed
      print("NOT REMOVING WHOLE POLYGON")
    }
    if (is.null(intList)) {
      #print("removing subsumed polygons")
      polytrim = rgeos::gUnion(polytrim, polytrim)
    } else {
      for (int in 1:length(intList)) {
        polytrim = rgeos::gDifference(polytrim, intList[[int]])
      }
      #print("removing subsumed polygons and differences from larger shape")
    }
  }
  return(polytrim)
}

