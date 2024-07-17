#' @import raster
#' @import MASS
#' @import spocc
#' @import rgeos
#' @import dplyr
#' @import sp
#' @import viridis
#' @import caret
#' @import sf
#' @import rebird
NULL
#' Flag Bad Polygon Overlaps 2
#'
#' This function takes two polygons representing subspecies ranges for Subspecies-A and Subspecies-B
#' and flags for removal intersections between the polygons. Called within polygonTrimmer()
#'
#' Checks for intersections in two ways:
#'
#' 1) Any polygon features that are completely subsumed within the other
#' subspecies range are flagged for removal, e.g., if Subspecies-A feature is completely subsumed by
#' Subspecies-B, then the Subspecies-A feature is flagged for removal in its entirety.
#' 2) Any intersections between polygons that do not comprise entire polygon features are flagged
#' for removal such that the intersection is flagged for removal from the larger polygon, but
#' not the smaller polygon. If polygons are the same size, intersection is removed from both.
#' E.g., Subspecies-A has a larger range than Subspecies-B. An intersection between them is removed
#' only from Subspecies-A.
#'
#' @param subsppPoly1 Polygon for first subspecies to compare
#' @param subsppPoly2 Polygon for second subspecies to compare
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
flagPolygonOverlap2 = function(subsppPoly1 = polA,
                               subsppPoly2 = polB,
                               crs = "+proj=longlat +ellps=WGS84") {
  ##function(subsppPoly1,subsppPoly2)
  ## this function checks for overlaps between polygons
  ## TODO: remove polygon if not touching another of same spp
  ## but also closer to polygon of other species
  #library(rgeos)
  #library(raster)
  ## there is a bug -- if one subspp range is entirely subsumed within another polygon,
  ## will delete that subspecies. no bueno
  badList_subsppA_features = c()
  badList_subsppB_features = c()
  overlapsToRemove_subsppA = c()
  overlapsToRemove_subsppB = c()
  ## NOT RIGHT BECAUSE NOT ADDING TO CORRECT SPOTS
  for (feature_subsppA in (1:length(subsppPoly1))) {
    ## get the features within the subspecies polygon
    for (feature_subsppB in (1:length(subsppPoly2))) {
      ## get the features within the subspecies polygon
      ## check areas
      ## Warning message:
      ## In RGEOSMiscFunc(spgeom, byid, "rgeos_area") :
      ##  Spatial object is not projected; GEOS expects planar coordinates
      totArea1 = rgeos::gArea(subsppPoly1) ## the whole area of the subspecies
      totArea2 = rgeos::gArea(subsppPoly2) ## the whole area of the subspecies
      area1 = rgeos::gArea(subsppPoly1[feature_subsppA,]) ## the area of the single feature
      area2 = rgeos::gArea(subsppPoly2[feature_subsppB,]) ## the area of the single feature
      # subsppPoly1$totalArea = totArea1
      # subsppPoly2$totalArea = totArea2
      # subsppPoly1[feature_subsppA,]$featureArea = area1
      # subsppPoly2[feature_subsppB]$featureArea = area2
      if (rgeos::gIntersects(subsppPoly1[feature_subsppA,], subsppPoly2[feature_subsppB,])) {
        #print("INTERSECTS")
        testInt = rgeos::gIntersection(subsppPoly1[feature_subsppA,], subsppPoly2[feature_subsppB,])
        #intersect = gIntersection(subsppPoly1[feature_subsppA,],subsppPoly2[feature_subsppB,],byid=T)
        #print(length(intersect))
        ## if they overlap
        #plot(subsppPoly1[feature_subsppA,],border="red",add=T)
        #plot(subsppPoly2[feature_subsppB,],border="cyan",add=T)
        testArea = rgeos::gArea(testInt)
        if (testArea > 0) {
          #print("OVERLAP AREA NOT ZERO")
          intersect = raster::intersect(subsppPoly1[feature_subsppA,], subsppPoly2[feature_subsppB,])
          #intersect = gIntersection(subsppPoly1[feature_subsppA,],subsppPoly2[feature_subsppB,],byid=T)
          areaInt = rgeos::gArea(intersect)
          ## if they overlap
          area1percent = areaInt / totArea1 ## the area of the feature as a percent of the area of the subspecies
          area2percent = areaInt / totArea2
          ## NEW IF STATEMENTS
          if (area1 > area2) {
            ## if A big B small -- area1 vs area2
            if (testArea < area2) {
              ## if B not subsumed, remove from A -- testArea
              #badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
              overlapsToRemove_subsppA = c(overlapsToRemove_subsppA, intersect)
            }
            else if (testArea >= area2) {
              ## else if B subsumed
              if (area2percent == 1) {
                ## if B = 100% of subspecies, remove from A -- area2percent
                #badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
                overlapsToRemove_subsppA = c(overlapsToRemove_subsppA, intersect)
              }
              else if (area2percent != 1) {
                ## else if B not 100% of subspecies, remove from B
                badList_subsppB_features = c(badList_subsppB_features, feature_subsppB)
              }
            }
          }
          else if (area1 < area2) {
            ## else if A small B big
            if (testArea < area1) {
              ## if A not subsumed, remove from B -- area1percent
              #badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
              overlapsToRemove_subsppB = c(overlapsToRemove_subsppB, intersect)
            }
            else if (testArea >= area1) {
              ## else if A subsumed
              if (area1percent == 1) {
                ## if A = 100% of subspecies, remove from B
                #badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
                overlapsToRemove_subsppB = c(overlapsToRemove_subsppB, intersect)
              }
              else if (area1percent != 1) {
                ## else if A not 100% of subspecies, remove from A
                badList_subsppA_features = c(badList_subsppA_features, feature_subsppA)
              }
            }
          }
          else if (area1 == area2) {
            ## else if A = B
            ## check total areas and remove larger
            if (totArea1 > totArea2) {
              ## if total A is greater than total B, remove A
              #badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
              overlapsToRemove_subsppA = c(overlapsToRemove_subsppA, intersect)
            }
            else if (totArea1 < totArea2) {
              ## if total B is greater than total A, remove B
              #badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
              overlapsToRemove_subsppB = c(overlapsToRemove_subsppB, intersect)
            }
            else if (totArea1 == totArea2) {
              ## if they are the exact same size
              if (testArea < area1 &&
                  testArea < area2) {
                ## if neither subsumed, remove arbitrary
                flip = sample(1:2, 1)
                if (flip == 1) {
                  ## flip coin remove A
                  #badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
                  overlapsToRemove_subsppA = c(overlapsToRemove_subsppA, intersect)
                }
                else if (flip == 2) {
                  ## flip coin remove B
                  #badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
                  overlapsToRemove_subsppB = c(overlapsToRemove_subsppB, intersect)
                }
              }
              else if (testArea == area1 &&
                       testArea == area2) {
                ## if both are subsumed, keep both
                print("EXACT MATCH, KEEPING BOTH")
              }
            }
          }
        }
      }
    }
    toReturn = list(
      subsppApoly_toremove = badList_subsppA_features,
      subsppBpoly_toremove = badList_subsppB_features,
      subsppA_intToRemove = overlapsToRemove_subsppA,
      subsppB_intToRemove = overlapsToRemove_subsppB
    )
    return(toReturn)
  }
}

