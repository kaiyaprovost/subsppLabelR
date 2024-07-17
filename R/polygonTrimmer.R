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
#' Polygon Trimmer
#'
#' This function takes a list of polygons as from densityMapToPolygons() or similar
#' and removes portions of the polygons as determined by flagPolygonOverlap() and
#' trimPolygonsByOverlap()
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
#' densityPolygons_trim = polygonTrimmer(polygonList=densityPolygons,
#'    namesList=c("sinuatus","fulvescens"))
#'
polygonTrimmer = function(polygonList, namesList) {
  newPolygonList = polygonList
  for (slotA in 1:length(namesList)) {
    for (slotB in 1:length(namesList)) {
      if (namesList[[slotA]] != "unknown" &&
          namesList[[slotB]] != "unknown" && slotA != slotB) {
        print(paste(slotA,slotB,sep=" "))
        #print(paste(namesList[[slotA]],"with",namesList[[slotB]],sep=" "))
        polA = newPolygonList[[slotA]]
        polB = newPolygonList[[slotB]]
        #print(class(polA))
        #print(class(polB))
        # plot(bg,col="grey",colNA="darkgrey")
        # plot(polA,add=T,border="cyan",lwd=7)
        # plot(polB,add=T,border="red",lwd=4)
        # #invisible(readline(prompt="Press [enter] to continue"))
        # if(!is.null(raster::intersect(polA,polB))){
        #   plot(raster::intersect(polA,polB),add=T,lwd=1,border="black")
        # }
        ## this throws the warnings
        polygonsToRemove = (flagPolygonOverlap2(polA, polB)) ######### CHANGED
        #print("CALLING NEW FUNCTION")
        subsppA_polygonsToRemove = polygonsToRemove$subsppApoly_toremove
        subsppB_polygonsToRemove = polygonsToRemove$subsppBpoly_toremove
        overlapToRemove_subsppA = polygonsToRemove$subsppA_intToRemove
        overlapToRemove_subsppB = polygonsToRemove$subsppB_intToRemove
        subsppA_densityPolygon_trim = trimPolygonsByOverlap(polygon = polA,
                                                            idList = subsppA_polygonsToRemove,
                                                            intList = overlapToRemove_subsppA)
        ## THIS IS FAILING
        subsppB_densityPolygon_trim = trimPolygonsByOverlap(polygon = polB,
                                                            idList = subsppB_polygonsToRemove,
                                                            intList = overlapToRemove_subsppB)
        ## TODO: fix this
        subsppA = namesList[[slotA]]
        subsppB = namesList[[slotB]]
        if(!(is.null(subsppA_densityPolygon_trim))){newPolygonList[[slotA]] = subsppA_densityPolygon_trim}
        if(!(is.null(subsppB_densityPolygon_trim))){newPolygonList[[slotB]] = subsppB_densityPolygon_trim}
        names(newPolygonList) = names(polygonList)
        #print(names(newPolygonList))
        #plot(bg,col="grey",colNA="darkgrey")
        #plot(subsppA_densityPolygon_trim,add=T,border="cyan",lwd=7)
        #plot(subsppB_densityPolygon_trim,add=T,border="red",lwd=4)
        #plot(intersect(polA,polB),add=T,lwd=1,border="black")
      }
    }
  }
  return(newPolygonList)
}
