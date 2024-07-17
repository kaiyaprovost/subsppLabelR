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
#' Check for Subspecies Contiguity
#'
#' This function takes multiple subspecies ranges as polygons and checks whether
#' 1) polygon features are adjacent to their own subspecies, and if not, 2)
#' whether polygon features are closer to their own subspecies than to another.
#' If a subspecies feature does not meet either criteria, it is eliminated.
#'
#' WARNING!!! Function is not working. Do not use.
#'
#' @param listOfPolygons List of subspecies range polygons.
#'
#' @export
#'
closestPolygonFunction = function(listOfPolygons) {
  ## TODO: is broken
  ##TODO: figure out how to do so that you check if polygons are between other polygons
  ## for now just remove if touching non-self boundary
  ## check each feature of each polygon
  newListPoly = listOfPolygons
  polyNames = names(newListPoly)
  for (subspp_i in 1:length(newListPoly)) {
    if (polyNames[subspp_i] != "unknown") {
      print(polyNames[subspp_i])
      wholePol = newListPoly[[subspp_i]]
      wholePol = sp::disaggregate(wholePol)
      if (length(wholePol) <= 1) {
        print("only one feature, skipping")
      }
      else {
        wholePol$area = rgeos::gArea(wholePol, byid = T)
        print(wholePol$area)
        ## sort features by size
        wholePol = wholePol[rev(order(wholePol$area)),]
        #print(wholePol)
        for (features_j in 1:length(wholePol)) {
          feat2check = wholePol[features_j,]
          if (rgeos::gArea(feat2check) < max(wholePol$area)) {
            ## if not the largest chunk, check if touching other chunks
            int = (rgeos::gIntersection(feat2check, wholePol[-features_j,]))
            if (is.null(int)) {
              ## diagonal is fine
              ## if not touching, check what polygon is shortest distance to
              polyListForCheck = newListPoly
              polyListForCheck[[subspp_i]] = wholePol[-features_j,]
              #plot(newListPoly[[subspp_i]])
              #plot(polyListForCheck[[subspp_i]])
              #polyListForCheck = polyListForCheck[names(polyListForCheck)!="unknown"]
              dists = sapply(1:length(polyListForCheck), function(i) {
                x = ((
                  rgeos::gDistance(feat2check, polyListForCheck[[i]])
                ))
                names(x) = names(polyListForCheck)[[i]]
                return(x)
              })
              dists[which(names(dists) == "unknown")] = 1e9
              print(paste(
                "Testing ",
                names(polyListForCheck)[subspp_i],
                " against ",
                paste(names(which(
                  dists == min(dists)
                )), collapse = ", "),
                sep = ""
              ))
              if (!(polyNames[subspp_i] %in% names(which(dists == min(
                dists
              ))))) {
                print("Closest is not same subspecies")
                if (min(dists) == 0) {
                  print("Touching other subspecies, removing")
                  removed = wholePol[-features_j,]
                  newListPoly[[subspp_i]] = removed
                }
              }
            }
          }
        }
      }
    }
  }
  return(newListPoly)
}

