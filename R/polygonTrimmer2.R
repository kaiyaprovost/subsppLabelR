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
#' Polygon Trimmer 2
#'
#' This function takes a list of polygons as from densityMapToPolygons() or similar
#' and removes portions of the polygons based off of their overlaps.
#'
#' @param polygonList A list of polygons to check
#' @param namesList The subspecies (or other) names associated with polygon list
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#'
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
#' densityPolygons = list(sinuatus=densPol_sin,fulvescens=densPol_ful)
#' densityPolygons_trim = polygonTrimmer2(polygonList=densityPolygons,
#'    namesList=c("sinuatus","fulvescens"))
#'
polygonTrimmer2 = function(polygonList, namesList, crs = "+proj=longlat +ellps=WGS84") {
  ## if the names list has changed this doesn't work?
  newPolygonList = polygonList

  if(length(namesList)!=length(polygonList)) {
    namesList = names(polygonList)
    print(names(polygonList))
  }

  print(paste(length(newPolygonList),length(namesList)))

  for (slotA in 1:length(namesList)) {
    for (slotB in 1:length(namesList)) {
      if (namesList[[slotA]] != "unknown" &&
          namesList[[slotB]] != "unknown" && slotA != slotB) {
        print(paste(slotA,namesList[[slotA]],slotB,namesList[[slotB]],sep=" "))

        try({
          polA = newPolygonList[[slotA]]
          polB = newPolygonList[[slotB]]

          ## check overlap between polA and polB
          overlapPol = sf::st_intersection(polA,polB)

          if(is.null(overlapPol)) {
            overlapArea = 0
          } else {
            overlapArea = sf::st_Area(overlapPol)
          }

          if(overlapArea != 0){
            if(class(overlapPol)=="SpatialCollections") {
              overlapPol = overlapPol@polyobj
            }
            ## check the overlap size relative to the other sizes
            areaPolA = sf::st_Area(polA)
            areaPolB = sf::st_Area(polB)

            ## check if it is smaller than one or the other or both

            if(overlapArea < areaPolA) {all_of_A = F} else {all_of_A = T}
            if(overlapArea < areaPolB) {all_of_B = F} else {all_of_B = T}

            if (all_of_A == F && all_of_B == T) {
              ## remove from A
              polA = sf::st_Difference(polA,overlapPol) ## order matters
            } else if (all_of_A == T && all_of_B == F) {
              ## remove from B
              polB = sf::st_Difference(polB,overlapPol) ## order matters

            } else if (all_of_A == F && all_of_B == F) {
              ## remove from both
              polA = sf::st_Difference(polA,overlapPol) ## order matters
              polB = sf::st_Difference(polB,overlapPol) ## order matters

            } ## we don't do anything if both are true


          }

          newPolygonList[[slotA]] = polA
          newPolygonList[[slotB]] = polB
        })

      }
    }
  }
  return(newPolygonList)
}
